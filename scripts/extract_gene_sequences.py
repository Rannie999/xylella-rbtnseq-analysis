#!/usr/bin/env python3
"""Map locusIds from top/bottom CSVs to gene coordinates and extract FASTA."""
import argparse, csv, re, sys
from pathlib import Path

ALIASES = {
 "locus": ("locusid", "locus_id", "geneid", "gene_id", "id"),
 "scaffold": ("scaffold", "scaffoldid", "scaffold_id", "seqid", "chr", "chrom", "chromosome"),
 "start": ("start", "gene_start", "begin"), "end": ("end", "gene_end", "stop"),
 "strand": ("strand",),
}

def norm(s): return s.strip().lower().replace(" ", "").replace("-", "_")

def reader(path):
    h = path.open(encoding="utf-8-sig", newline="")
    sample = h.read(16384); h.seek(0)
    try: dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error: dialect = csv.excel
    return h, csv.DictReader(h, dialect=dialect)

def column(fields, kind, required=True):
    lookup = {norm(x): x for x in fields or []}
    value = next((lookup[x] for x in ALIASES[kind] if x in lookup), None)
    if required and value is None:
        raise ValueError(f"No {kind} column; available: {', '.join(fields or [])}")
    return value

def locus_ids(path):
    h, rows = reader(path)
    with h:
        col = column(rows.fieldnames, "locus"); result=[]; seen=set()
        for row in rows:
            value=(row.get(col) or "").strip()
            if value and value not in seen: result.append(value); seen.add(value)
    return result

def locations(path):
    h, rows = reader(path)
    with h:
        fields=rows.fieldnames
        lc=column(fields,"locus"); sc=column(fields,"scaffold")
        st=column(fields,"start"); en=column(fields,"end")
        sr=column(fields,"strand",False); result={}
        for number,row in enumerate(rows,2):
            locus=(row.get(lc) or "").strip()
            if not locus: continue
            try: start=int(row[st]); end=int(row[en])
            except (ValueError,TypeError): raise ValueError(f"Bad coordinates at {path}:{number}")
            loc=((row.get(sc) or "").strip(),start,end,(row.get(sr) or "+").strip() if sr else "+")
            if locus in result and result[locus] != loc: raise ValueError(f"Conflicting rows for {locus}")
            result[locus]=loc
    return result

def fasta(path):
    result={}; key=None; parts=[]
    def save():
        if key is not None:
            if key in result: raise ValueError(f"Duplicate FASTA ID: {key}")
            result[key]="".join(parts).upper()
    with path.open(encoding="utf-8") as h:
        for n,raw in enumerate(h,1):
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                save(); key=line[1:].split()[0]; parts=[]
            elif key is None: raise ValueError(f"Sequence before FASTA header at line {n}")
            else: parts.append("".join(line.split()))
    save(); return result

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn","TGCAYRKMVHDBNtgcayrkmvhdbn"))[::-1]

def scaffold_lookup(genome):
    """Map both exact FASTA IDs and IDs without a final version suffix."""
    lookup={name:name for name in genome}
    versionless={}
    for name in genome:
        base=re.sub(r"\.\d+$", "", name)
        versionless.setdefault(base, []).append(name)
    return lookup, versionless

def resolve_scaffold(name, exact, versionless):
    if name in exact:
        return exact[name]
    matches=versionless.get(name, [])
    if len(matches)==1:
        return matches[0]
    if len(matches)>1:
        raise ValueError(
            f"Scaffold {name!r} matches multiple FASTA records: {', '.join(matches)}"
        )
    raise ValueError(f"Scaffold {name!r} absent from genome FASTA")

def write(label, ids, locs, genome, coord, outdir):
    missing=[]; count=0; exact,versionless=scaffold_lookup(genome)
    with (outdir/f"{label}_gene_sequences.fasta").open("w") as out:
        for locus in ids:
            if locus not in locs: missing.append(locus); continue
            scaffold,start,end,strand=locs[locus]
            fasta_scaffold=resolve_scaffold(scaffold,exact,versionless)
            begin=start-1 if coord=="one-based-inclusive" else start
            if begin < 0 or end > len(genome[fasta_scaffold]): raise ValueError(f"Coordinates outside {fasta_scaffold} for {locus}")
            seq=genome[fasta_scaffold][begin:end]
            if strand == "-": seq=revcomp(seq)
            out.write(f">{locus} {fasta_scaffold}:{start}-{end}({strand})\n")
            out.write("\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n"); count+=1
    (outdir/f"{label}_missing_locusIds.txt").write_text("".join(x+"\n" for x in missing))
    print(f"{label}: wrote {count}/{len(ids)} sequences; {len(missing)} locusIds missing")
    return not missing

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top",type=Path,default=Path("top_gene.csv"))
    p.add_argument("--bottom",type=Path,default=Path("bottom_genes.csv"))
    p.add_argument("--genes",type=Path,default=Path("exp4/time0_B/genes"))
    p.add_argument("--genome",type=Path,required=True,help="scaffold/genome FASTA")
    p.add_argument("--output-dir",type=Path,default=Path("extracted_sequences"))
    p.add_argument("--coordinates",choices=("one-based-inclusive","zero-based-half-open"),default="one-based-inclusive")
    a=p.parse_args()
    try:
        for path in (a.top,a.bottom,a.genes,a.genome):
            if not path.is_file(): raise ValueError(f"File not found: {path}")
        locs=locations(a.genes); genome=fasta(a.genome); a.output_dir.mkdir(parents=True,exist_ok=True)
        ok1=write("top",locus_ids(a.top),locs,genome,a.coordinates,a.output_dir)
        ok2=write("bottom",locus_ids(a.bottom),locs,genome,a.coordinates,a.output_dir)
        return 0 if ok1 and ok2 else 1
    except (OSError,ValueError,csv.Error) as e: print(f"Error: {e}",file=sys.stderr); return 2

if __name__ == "__main__": raise SystemExit(main())
