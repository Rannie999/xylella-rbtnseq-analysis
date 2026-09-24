#!/usr/bin/env python3
"""Summarize barcode coverage and fitness results for one or more FEBA runs."""
from __future__ import annotations

import argparse, csv, math, re
from pathlib import Path

ID_COLUMNS = {"locusId", "sysName", "desc"}

def args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run_dir", required=True, type=Path,
                   help="one completed run directory, or a runs/ directory containing completed runs")
    return p.parse_args()

def tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def number(value):
    try:
        value = float(value or "")
        return value if math.isfinite(value) else None
    except ValueError:
        return None

def mean(values): return sum(values) / len(values) if values else None
def pct(n, d): return f"{100 * n / d:.2f}%" if d else "NA"
def token(value):
    match = re.search(r"S\d+", value)
    return match.group(0) if match else None

def quality(feba):
    path = feba / "fit_quality.tab"
    return tsv(path) if path.is_file() else []

def barcode_metrics(run):
    """Aggregate the run's generated barcode-QC results.

    This deliberately uses barcode_qc.py outputs so coverage, insertion, and
    retention values have the same definitions in the per-sample and multi-run
    reports.
    """
    qc = run / "analyses" / "barcode_qc"
    sample_rows = list(csv.DictReader((qc / "sample_qc.csv").open(encoding="utf-8-sig", newline="")))
    gene_rows = list(csv.DictReader((qc / "gene_qc.csv").open(encoding="utf-8-sig", newline="")))
    if not sample_rows or not gene_rows:
        raise ValueError(f"Missing or empty barcode-QC results in {qc}")
    if len(sample_rows) != len(gene_rows):
        raise ValueError(f"Sample and gene QC row counts differ in {qc}")

    # Barcode-QC uses names beginning with T0 for the Time0 reference samples.
    selected = [r for r in sample_rows if not r.get("sample", "").startswith("T0")]
    retention_mean = mean([number(r.get("t0_barcode_retention_pct")) or 0 for r in selected])
    return {
        "Samples": len(sample_rows),
        "Mean detected barcodes per sample": round(mean([number(r.get("detected_barcodes")) or 0 for r in sample_rows]) or 0, 2),
        "Mean detected genes per sample": round(mean([number(r.get("detected_genes")) or 0 for r in gene_rows]) or 0, 2),
        "Mean reads per sample": round(mean([number(r.get("total_reads")) or 0 for r in sample_rows]) or 0, 2),
        "Mean insertions per gene": round(mean([number(r.get("mean_insertions_per_gene")) or 0 for r in gene_rows]) or 0, 2),
        "Mean Time0 barcode retention": f"{retention_mean:.2f}%" if retention_mean is not None else "NA",
    }

def fitness_metrics(run, threshold):
    feba = run / "feba"; fit, fit_t = tsv(feba / "fit_logratios.tab"), tsv(feba / "fit_t.tab")
    if not fit or not fit_t: raise ValueError(f"Missing or empty fitness tables in {feba}")
    samples = [c for c in fit[0] if c not in ID_COLUMNS and "time0" not in c.lower()]
    t_rows = {r.get("locusId", ""): r for r in fit_t}
    values, t_values, genes, passing_genes = [], [], set(), set()
    above_fit = above_t = 0
    for row in fit:
        locus, t_row, represented = row.get("locusId", ""), t_rows.get(row.get("locusId", ""), {}), False
        for sample in samples:
            value, t_value = number(row.get(sample)), number(t_row.get(sample))
            if value is not None:
                represented = True; values.append(value); above_fit += abs(value) > 1
            if t_value is not None:
                t_values.append(t_value)
                if abs(t_value) > threshold: above_t += 1; passing_genes.add(locus)
        if represented: genes.add(locus)
    possible = len(genes) * len(samples)
    q = quality(feba); experimental = [r for r in q if r.get("short", "").strip().lower() != "time0"]
    passed = [r for r in experimental if r.get("u", "").upper() == "TRUE"]
    return {
        "Samples analyzed": len(samples), "Genes represented": len(genes),
        "Possible gene-sample values": possible, "Non-missing fitness values": len(values),
        "Gene-sample values with |fitness| > 1": f"{above_fit:,} ({pct(above_fit, possible)})",
        f"Values passing |t-like| > {threshold:g}": f"{above_t:,} ({pct(above_t, possible)})",
        "Genes with at least one passing value": len(passing_genes),
        "Overall mean absolute fitness": round(mean([abs(x) for x in values]) or 0, 4),
        "Overall mean absolute t-like": round(mean([abs(x) for x in t_values]) or 0, 4),
        "FEBA quality-passing experimental samples": f"{len(passed)}/{len(experimental)}",
    }

def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=fields); out.writeheader(); out.writerows(rows)

def write_md(path, coverage, fitness):
    def cell(value): return str(value).replace("|", "\\|")
    columns = list(coverage[0])
    text = ["# FEBA multi-run summary", "", "## Barcode and gene coverage", "",
            "| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    text += ["| " + " | ".join(cell(row[c]) for c in columns) + " |" for row in coverage]
    columns = list(fitness[0])
    text += ["", "## Fitness results", "", "| " + " | ".join(columns) + " |",
             "| " + " | ".join(["---"] * len(columns)) + " |"]
    text += ["| " + " | ".join(cell(row[c]) for c in columns) + " |" for row in fitness]
    path.write_text("\n".join(text) + "\n", encoding="utf-8")

def main():
    a = args(); selected = a.run_dir.resolve()
    if (selected / "feba").is_dir():
        runs = [selected]
    elif selected.is_dir():
        runs = sorted(path for path in selected.iterdir() if path.is_dir() and (path / "feba").is_dir())
    else:
        raise SystemExit(f"Run directory not found: {selected}")
    if not runs:
        raise SystemExit(f"No completed run directories containing feba/ found in: {selected}")
    labels = [path.name for path in runs]
    for run in runs:
        for name in ("all.poolcount", "fit_logratios.tab", "fit_t.tab"):
            if not (run / "feba" / name).is_file(): raise SystemExit(f"Missing: {run / 'feba' / name}")
    if len(runs) == 1: output = runs[0] / "summaries"
    else:
        output = selected / "multi_run_summary"
    output.mkdir(parents=True, exist_ok=True)
    coverage = [{"Experiment": label, **barcode_metrics(run)} for run, label in zip(runs, labels)]
    data = {label: fitness_metrics(run, 2.5) for run, label in zip(runs, labels)}
    fitness = [{"Experiment": label, **data[label]} for label in labels]
    write_csv(output / "barcode_gene_summary.csv", coverage, list(coverage[0]))
    write_csv(output / "fitness_summary.csv", fitness, list(fitness[0]))
    write_md(output / "analysis_summary.md", coverage, fitness)
    print(f"Wrote barcode/gene summary: {output / 'barcode_gene_summary.csv'}")
    print(f"Wrote fitness summary: {output / 'fitness_summary.csv'}")
    print(f"Wrote report: {output / 'analysis_summary.md'}")

if __name__ == "__main__": main()
