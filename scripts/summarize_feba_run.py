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
    feba = run / "feba"; rows = tsv(feba / "all.poolcount")
    if not rows: raise ValueError(f"No barcode rows in {feba / 'all.poolcount'}")
    fields = {"barcode", "rcbarcode", "scaffold", "strand", "pos", "locusId", "f"}
    samples = [col for col in rows[0] if col not in fields]
    detected = [r for r in rows if any((number(r.get(c)) or 0) > 0 for c in samples)]
    genic = [r for r in detected if r.get("locusId", "").strip()]
    genes = {r["locusId"].strip() for r in genic}
    time0 = {token(r.get("name", "")) for r in quality(feba)
             if r.get("short", "").strip().lower() == "time0"}
    time0.discard(None)
    time0_barcodes = {r["barcode"] for r in detected if any(
        (number(r.get(c)) or 0) > 0 and token(c) in time0 for c in samples)}
    barcode_per, gene_per, reads_per, retention = [], [], [], []
    for sample in samples:
        present = [r for r in detected if (number(r.get(sample)) or 0) > 0]
        barcode_per.append(len(present))
        gene_per.append(len({r["locusId"].strip() for r in present if r.get("locusId", "").strip()}))
        reads_per.append(sum(number(r.get(sample)) or 0 for r in present))
        if token(sample) not in time0 and present and time0_barcodes:
            retention.append(sum(r["barcode"] in time0_barcodes for r in present) / len(present))
    return {
        "Samples": len(samples), "Detected barcodes": len(detected), "Detected genes": len(genes),
        "Genic barcodes": len(genic), "Intergenic barcodes": len(detected) - len(genic),
        "Mean barcodes per sample": round(mean(barcode_per) or 0, 2),
        "Mean genes per sample": round(mean(gene_per) or 0, 2),
        "Mean reads per sample": round(mean(reads_per) or 0, 2),
        "Mean insertions per gene": round(len(genic) / len(genes), 2) if genes else "NA",
        "Mean Time0 barcode retention": pct(round(10000 * (mean(retention) or 0)), 10000) if retention else "NA",
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
