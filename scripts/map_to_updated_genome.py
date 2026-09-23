#!/usr/bin/env python3
"""BLAST old gene sequences onto a new genome and recover loci/DNA/proteins."""
import argparse, csv, hashlib, re, shutil, subprocess, sys
from pathlib import Path
from urllib.parse import unquote

def fasta(path):
    records={}; key=None; header=None; parts=[]
    def save():
        if key: records[key]=(header,"".join(parts).upper())
    with path.open(encoding="utf-8") as f:
        for raw in f:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                save(); header=line[1:]; key=header.split()[0]; parts=[]
            elif key: parts.append("".join(line.split()))
            else: raise ValueError(f"{path} is not FASTA-formatted")
    save(); return records

def write_fasta_record(out, header, seq):
    out.write(">"+header+"\n")
    out.write("\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n")

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTRYMKBDHVN","TGCAYRKMVHDBN"))[::-1]

def attrs(text):
    result={}
    for item in text.split(";"):
        if "=" in item:
            k,v=item.split("=",1); result[k]=unquote(v)
    return result

def gff(path):
    features=[]
    with path.open(encoding="utf-8") as f:
        for n,line in enumerate(f,1):
            if not line.strip() or line.startswith("#"): continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)!=9: raise ValueError(f"Invalid GFF3 row at {path}:{n}")
            features.append({"seqid":fields[0],"type":fields[2],"start":int(fields[3]),
                "end":int(fields[4]),"strand":fields[6],"attrs":attrs(fields[8])})
    genes=[x for x in features if x["type"].lower() in ("gene","pseudogene")]
    return features,genes

def base_id(value): return re.sub(r"\.\d+$","",value)

def seqid_map(genome):
    bybase={}
    for x in genome: bybase.setdefault(base_id(x),[]).append(x)
    return bybase

def resolve_seqid(value, genome, bybase):
    if value in genome: return value
    choices=bybase.get(base_id(value),[])
    if len(choices)==1: return choices[0]
    raise ValueError(f"Cannot uniquely match scaffold {value!r} to genome FASTA")

def run_blast(query, genome, out, threads):
    for tool in ("makeblastdb","blastn"):
        if not shutil.which(tool): raise ValueError(f"{tool} not found; install NCBI BLAST+")
    db=out.parent/"blastdb"/"updated_genome"
    db.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(["makeblastdb","-in",str(genome),"-dbtype","nucl","-out",str(db)],check=True)
    fields="qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore"
    subprocess.run(["blastn","-query",str(query),"-db",str(db),"-out",str(out),
        "-outfmt","6 "+fields,"-max_target_seqs","20","-num_threads",str(threads)],check=True)

def best_hits(path, min_identity, min_coverage):
    hits={}
    with path.open() as f:
        for row in csv.reader(f,delimiter="\t"):
            q,s=row[0],row[1]; ident=float(row[2]); length=int(row[3]); qlen=int(row[4])
            hit={"query":q,"scaffold":s,"identity":ident,"length":length,"qlen":qlen,
                "qstart":int(row[5]),"qend":int(row[6]),"sstart":int(row[7]),"send":int(row[8]),
                "evalue":float(row[9]),"bitscore":float(row[10])}
            hit["coverage"]=100.0*length/qlen
            if ident < min_identity or hit["coverage"] < min_coverage: continue
            old=hits.get(q)
            if old is None or (hit["bitscore"],hit["coverage"],hit["identity"]) > (old["bitscore"],old["coverage"],old["identity"]): hits[q]=hit
    return hits

def overlap(a,b,c,d): return max(0,min(b,d)-max(a,c)+1)

def gene_for(hit, genes):
    lo,hi=sorted((hit["sstart"],hit["send"]))
    candidates=[]
    for gene in genes:
        if base_id(gene["seqid"])==base_id(hit["scaffold"]):
            ov=overlap(lo,hi,gene["start"],gene["end"])
            if ov: candidates.append((ov,gene))
    return max(candidates,key=lambda x:x[0])[1] if candidates else None

def feature_ids(gene, features):
    """Return identifiers only from this gene and its GFF3 descendants.

    Coordinate overlap is not sufficient: neighboring/nested genes can overlap and
    would otherwise lend each other their CDS protein_id.
    """
    if not gene: return []
    gene_id=gene["attrs"].get("ID","")
    descendants=[]
    if gene_id:
        children={}
        for feature in features:
            for parent in feature["attrs"].get("Parent","").split(","):
                if parent: children.setdefault(parent,[]).append(feature)
        queue=[gene_id]; seen={gene_id}
        while queue:
            parent=queue.pop(0)
            for child in children.get(parent,[]):
                child_id=child["attrs"].get("ID","")
                descendants.append(child)
                if child_id and child_id not in seen:
                    seen.add(child_id); queue.append(child_id)
    related=[gene]+descendants
    # Explicit protein_id values are the safest match and must be tried first.
    values=[]
    for feature in related:
        values.extend(feature["attrs"].get("protein_id","").split(","))
    for feature in related:
        for key in ("ID","Name","locus_tag","gene"):
            values.extend(feature["attrs"].get(key,"").split(","))
    result=[]; seen_values=set()
    for value in values:
        if value and value not in seen_values:
            result.append(value); seen_values.add(value)
    return result

def protein_for(ids, proteins):
    for value in ids:
        options=(value,base_id(value),value.removeprefix("gene:"),value.removeprefix("transcript:"),value.removeprefix("rna:"))
        for option in options:
            if option in proteins: return option
    for key,(header,seq) in proteins.items():
        if any(re.search(r"(?<![\w.-])"+re.escape(x)+r"(?![\w.-])",header) for x in ids): return key
    return None

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--query",type=Path,required=True,help="old extracted gene FASTA")
    p.add_argument("--new-genome",type=Path,required=True)
    p.add_argument("--gff3",type=Path,required=True)
    p.add_argument("--proteins",type=Path,required=True,help="FASTA amino acids; extension may be .txt")
    p.add_argument("--output-dir",type=Path,default=Path("updated_genome_results"))
    p.add_argument("--blast-results",type=Path,help="reuse existing BLAST tabular output")
    p.add_argument("--min-identity",type=float,default=70.0); p.add_argument("--min-coverage",type=float,default=70.0)
    p.add_argument("--threads",type=int,default=4); a=p.parse_args()
    try:
        a.output_dir.mkdir(parents=True,exist_ok=True)
        blastfile=a.blast_results or a.output_dir/"blast_hits.tsv"
        if not a.blast_results: run_blast(a.query,a.new_genome,blastfile,a.threads)
        queries=fasta(a.query); genome=fasta(a.new_genome); proteins=fasta(a.proteins)
        features,genes=gff(a.gff3); hits=best_hits(blastfile,a.min_identity,a.min_coverage); bybase=seqid_map(genome)
        table=a.output_dir/"new_locus_info.csv"
        protein_queries={}; sequence_assignments={}
        with table.open("w",newline="") as tf, (a.output_dir/"new_locus_sequences.fasta").open("w") as nf, (a.output_dir/"new_proteins.faa").open("w") as pf:
            cols=["query_locusId","blast_scaffold","blast_start","blast_end","blast_strand","identity","query_coverage","evalue","new_gene_id","new_gene_name","gene_scaffold","gene_start","gene_end","gene_strand","protein_id","status"]
            w=csv.DictWriter(tf,fieldnames=cols); w.writeheader()
            for query in queries:
                row={x:"" for x in cols}; row["query_locusId"]=query; hit=hits.get(query)
                if not hit: row["status"]="no_passing_blast_hit"; w.writerow(row); continue
                lo,hi=sorted((hit["sstart"],hit["send"])); strand="+" if hit["sstart"]<=hit["send"] else "-"; gene=gene_for(hit,genes)
                row.update(blast_scaffold=hit["scaffold"],blast_start=lo,blast_end=hi,blast_strand=strand,identity=f'{hit["identity"]:.3f}',query_coverage=f'{hit["coverage"]:.3f}',evalue=hit["evalue"])
                sid=resolve_seqid(hit["scaffold"],genome,bybase)
                if gene:
                    ids=feature_ids(gene,features); gid=gene["attrs"].get("ID",""); pname=protein_for(ids,proteins)
                    row.update(new_gene_id=gid,new_gene_name=gene["attrs"].get("Name",gene["attrs"].get("locus_tag","")),gene_scaffold=gene["seqid"],gene_start=gene["start"],gene_end=gene["end"],gene_strand=gene["strand"],protein_id=pname or "",status="matched_gene" if pname else "matched_gene_no_protein")
                    sid=resolve_seqid(gene["seqid"],genome,bybase); seq=genome[sid][1][gene["start"]-1:gene["end"]]
                    if gene["strand"]=="-": seq=revcomp(seq)
                    write_fasta_record(nf,f"{query} new_gene={gid} {sid}:{gene['start']}-{gene['end']}({gene['strand']})",seq)
                    if pname:
                        write_fasta_record(pf,f"{query} new_gene={gid} protein={pname}",proteins[pname][1])
                        protein_queries.setdefault(pname,[]).append(query)
                        sequence_assignments.setdefault(proteins[pname][1],[]).append((pname,query))
                else:
                    seq=genome[sid][1][lo-1:hi]; seq=revcomp(seq) if strand=="-" else seq
                    row["status"]="blast_hit_no_overlapping_gene"; write_fasta_record(nf,f"{query} blast_hit {sid}:{lo}-{hi}({strand})",seq)
                w.writerow(row)
        duplicate_report=a.output_dir/"duplicate_protein_assignments.tsv"
        with duplicate_report.open("w",newline="") as report:
            writer=csv.writer(report,delimiter="\t")
            writer.writerow(["duplicate_type","protein_ids","query_locusIds","sequence_sha256","note"])
            for protein_id,query_ids in protein_queries.items():
                if len(query_ids)>1:
                    writer.writerow(["same_protein_id",protein_id,",".join(query_ids),"","Review: queries mapped to the same updated protein/gene"])
            for sequence,assignments in sequence_assignments.items():
                if len(assignments)>1 and len({query for _,query in assignments})>1:
                    writer.writerow(["identical_amino_acid_sequence",",".join(sorted({protein for protein,_ in assignments})),
                        ",".join(query for _,query in assignments),hashlib.sha256(sequence.encode()).hexdigest(),
                        "May be valid if nucleotide differences are synonymous; review gene mappings"])
        print(f"Wrote results to {a.output_dir}"); return 0
    except (OSError,ValueError,subprocess.CalledProcessError,csv.Error) as e: print(f"Error: {e}",file=sys.stderr); return 2
if __name__=="__main__": raise SystemExit(main())
