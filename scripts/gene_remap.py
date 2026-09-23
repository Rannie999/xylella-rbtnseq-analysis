#!/usr/bin/env python3
"""One-command old-locus extraction and mapping to an updated genome."""
import argparse, csv, itertools, re, sys
from pathlib import Path
from urllib.parse import quote

import extract_gene_sequences as old
import map_to_updated_genome as new

FASTA_EXTENSIONS={".fa",".fasta",".fna",".fas"}

def base(value): return re.sub(r"\.\d+$","",value)

def clean_locus_id(value):
    return value[5:] if value.startswith("gene-") else value

def protein_name_from_faa_header(header):
    """Extract the NCBI [protein=...] value, or fall back to header text."""
    match=re.search(r"\[protein=([^\]]+)\]",header,flags=re.IGNORECASE)
    if match: return match.group(1).strip()
    parts=header.split(maxsplit=1)
    return parts[1].strip() if len(parts)>1 else ""

def one(paths, description):
    paths=list(paths)
    if len(paths)!=1:
        raise ValueError(f"Expected one {description}, found {len(paths)}: "+", ".join(str(x) for x in paths))
    return paths[0]

def discover_old_inputs(input_dir):
    positive=next((input_dir/x for x in ("top_positive.csv","top_positve.csv") if (input_dir/x).is_file()),None)
    if positive is None: raise ValueError(f"Missing top_positive.csv (or top_positve.csv) in {input_dir}")
    negative=input_dir/"top_negative.csv"
    for path in (positive,negative):
        if not path.is_file(): raise ValueError(f"Missing required file: {path}")
    return positive,negative

def discover_new_inputs(directory, require_proteins=True):
    gff3=one(itertools.chain(directory.glob("*.gff3"),directory.glob("*.gff")),"GFF3 file")
    preferred=("protein_csd.txt","protein_cds.txt","protein_cds.faa","proteins.faa","protein.faa")
    proteins=next((directory/x for x in preferred if (directory/x).is_file()),None)
    if proteins is None:
        candidates=[x for x in directory.iterdir() if x.is_file() and
                    any(word in x.name.lower() for word in ("protein","peptide","pep"))]
        if candidates: proteins=one(candidates,"protein FASTA")
        elif require_proteins: raise ValueError(f"No protein FASTA found in {directory}")
    fastas=[x for x in directory.iterdir() if x.is_file() and x.suffix.lower() in FASTA_EXTENSIONS
            and not any(word in x.name.lower() for word in ("protein","peptide","pep","cds"))]
    if not fastas: raise ValueError(f"No genomic .fa/.fasta/.fna file found in {directory}")
    if len(fastas)==1: genome=fastas[0]
    else:
        features,_=new.gff(gff3); wanted={base(x["seqid"]) for x in features}
        ranked=sorted(((len(wanted & {base(x) for x in new.fasta(path)}),path) for path in fastas),reverse=True,key=lambda x:x[0])
        if ranked[0][0]==0 or ranked[0][0]==ranked[1][0]:
            raise ValueError("Cannot uniquely identify the updated genomic FASTA from its GFF3 scaffold names")
        genome=ranked[0][1]
    return genome,gff3,proteins

def add_temecula1_annotations(primary_table, temecula_table, temeculaL_proteins_file):
    """Insert Temecula1 locus and TemeculaL FAA header as columns 3 and 4."""
    proteins=new.fasta(temeculaL_proteins_file)
    with temecula_table.open(encoding="utf-8-sig",newline="") as source:
        temecula={row["query_locusId"]:row for row in csv.DictReader(source)}
    temporary=primary_table.with_suffix(".temecula1.tmp")
    with primary_table.open(encoding="utf-8-sig",newline="") as source, temporary.open("w",encoding="utf-8",newline="") as target:
        rows=csv.DictReader(source)
        original=[x for x in (rows.fieldnames or []) if x not in
                  ("temecula1_locusId","Temecula1_locusId","temecula1_faa_header",
                   "Temecula1_faa_header","TemeculaL_faa_header","faa_header_description")]
        fields=original[:2]+["Temecula1_locusId","TemeculaL_faa_header"]+original[2:]
        writer=csv.DictWriter(target,fieldnames=fields); writer.writeheader()
        for row in rows:
            row.pop("faa_header_description",None)
            row.pop("Temecula1_faa_header",None)
            row.pop("temecula1_faa_header",None)
            match=temecula.get(row["query_locusId"],{})
            protein_id=(row.get("protein_id") or "").strip()
            protein=proteins.get(protein_id)
            if protein is None and protein_id:
                protein=next((record for key,record in proteins.items() if base(key)==base(protein_id)),None)
            row["Temecula1_locusId"]=clean_locus_id(match.get("new_gene_id",""))
            row["TemeculaL_faa_header"]=protein_name_from_faa_header(protein[0]) if protein else ""
            writer.writerow(row)
    temporary.replace(primary_table)

def add_comparison_locus(primary_table, comparison_table, column_name, after_column=None):
    with comparison_table.open(encoding="utf-8-sig",newline="") as source:
        loci={row["query_locusId"]:clean_locus_id(row.get("new_gene_id","")) for row in csv.DictReader(source)}
    temporary=primary_table.with_suffix(f".{column_name}.tmp")
    with primary_table.open(encoding="utf-8-sig",newline="") as source, temporary.open("w",encoding="utf-8",newline="") as target:
        rows=csv.DictReader(source); fields=[x for x in (rows.fieldnames or []) if x!=column_name]
        if after_column in fields: fields.insert(fields.index(after_column)+1,column_name)
        else: fields.append(column_name)
        writer=csv.DictWriter(target,fieldnames=fields); writer.writeheader()
        for row in rows:
            row[column_name]=loci.get(row["query_locusId"],""); writer.writerow(row)
    temporary.replace(primary_table)

def run_mapping(query, genome, gff3, proteins, output, args):
    saved=sys.argv
    sys.argv=["map_to_updated_genome.py","--query",str(query),"--new-genome",str(genome),
        "--gff3",str(gff3),"--proteins",str(proteins),"--output-dir",str(output),
        "--threads",str(getattr(args,"threads",4)),"--min-identity",str(getattr(args,"min_identity",70.0)),"--min-coverage",str(getattr(args,"min_coverage",70.0))]
    try: return new.main()
    finally: sys.argv=saved

def gene_annotations(genes_file):
    handle,rows=old.reader(genes_file)
    with handle:
        fields=rows.fieldnames or []
        locus_col=old.column(fields,"locus")
        lookup={old.norm(x):x for x in fields}
        desc_col=next((lookup[x] for x in ("description","desc","annotation","product","gene_description") if x in lookup),None)
        category_aliases=("functional_category","functionalcategory","function_category",
                          "functioncategory","category","cog_category","cogcategory","cog")
        category_col=next((lookup[x] for x in category_aliases if x in lookup),None)
        descriptions={}; categories={}
        for row in rows:
            locus=(row.get(locus_col) or "").strip()
            if locus:
                descriptions[locus]=(row.get(desc_col) or "").strip() if desc_col else ""
                category=(row.get(category_col) or "").strip() if category_col else ""
                categories[locus]=category or "Uncategorized"
        return descriptions,categories

def enrich_locus_table(table, positive_ids, negative_ids, descriptions, categories, proteins_file):
    positive=set(positive_ids); negative=set(negative_ids)
    proteins=new.fasta(proteins_file)
    temporary=table.with_suffix(".tmp")
    with table.open(encoding="utf-8-sig",newline="") as source, temporary.open("w",encoding="utf-8",newline="") as target:
        rows=csv.DictReader(source)
        added=("description","functional_category","fitness_change","paperblast_url")
        excluded=set(added) | {"protein_description","faa_header_description"}
        fields=["query_locusId",*added]+[x for x in (rows.fieldnames or []) if x!="query_locusId" and x not in excluded]
        writer=csv.DictWriter(target,fieldnames=fields); writer.writeheader()
        for row in rows:
            row.pop("protein_description",None)
            locus=row["query_locusId"]
            if locus in positive and locus in negative: change="positive_and_negative"
            elif locus in positive: change="positive"
            elif locus in negative: change="negative"
            else: change="unknown"
            protein_id=(row.get("protein_id") or "").strip()
            protein=proteins.get(protein_id)
            if protein is None and protein_id:
                protein=next((record for key,record in proteins.items() if base(key)==base(protein_id)),None)
            paperblast_url=""
            if protein:
                header,sequence=protein
                paperblast_url="https://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query="+quote(sequence,safe="")
            row["description"]=descriptions.get(locus,"")
            row["functional_category"]=categories.get(locus,"Uncategorized")
            row["fitness_change"]=change
            row["paperblast_url"]=paperblast_url
            writer.writerow(row)
    temporary.replace(table)

def remove_if_empty(path):
    """Delete an empty text output, including a header-only TSV report."""
    if not path.is_file(): return
    lines=[line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) <= 1: path.unlink()

def column_prefix(name):
    """Make a stable CSV-column prefix from an updated-genome directory name."""
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_") or "genome"

def protein_header_lookup(proteins_file):
    """Return protein IDs mapped to their full FAA headers."""
    proteins=new.fasta(proteins_file)
    headers={key: record[0] for key,record in proteins.items()}
    by_base={base(key): header for key,header in headers.items()}
    return headers,by_base

def write_summary_table(output, positive_ids, negative_ids, descriptions, mappings):
    """Combine locus assignments from every updated genome into one review table."""
    query_ids=list(dict.fromkeys([*positive_ids,*negative_ids]))
    positive=set(positive_ids); negative=set(negative_ids)
    fields=["query_locusId","description","fitness_change"]
    for mapping in mappings:
        fields.extend([f"{mapping['prefix']}_locusId",f"{mapping['prefix']}_faa_header"])
    fields.append("paperblast_url")

    for mapping in mappings:
        table=mapping["directory"] / "new_locus_info.csv"
        with table.open(encoding="utf-8-sig",newline="") as handle:
            mapping["rows"]={row["query_locusId"]:row for row in csv.DictReader(handle)}

    table=output/"new_locus_info.csv"
    with table.open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader()
        for locus in query_ids:
            if locus in positive and locus in negative: change="positive_and_negative"
            elif locus in positive: change="positive"
            elif locus in negative: change="negative"
            else: change="unknown"
            row={"query_locusId":locus,"description":descriptions.get(locus,""),
                 "fitness_change":change,"paperblast_url":""}
            for mapping in mappings:
                result=mapping["rows"].get(locus,{})
                protein_id=(result.get("protein_id") or "").strip()
                header=mapping["protein_headers"].get(protein_id)
                if header is None and protein_id:
                    header=mapping["protein_headers_by_base"].get(base(protein_id),"")
                prefix=mapping["prefix"]
                row[f"{prefix}_locusId"]=clean_locus_id(result.get("new_gene_id", ""))
                row[f"{prefix}_faa_header"]=header or ""
                if not row["paperblast_url"]:
                    row["paperblast_url"]=result.get("paperblast_url","")
            writer.writerow(row)
    print(f"Wrote combined locus summary: {table}")

def discover_genomes(genome_dir):
    """Find the old genome at the root and one updated genome in each subdirectory."""
    root_fastas=[path for path in genome_dir.iterdir() if path.is_file()
                 and path.suffix.lower() in FASTA_EXTENSIONS
                 and not any(word in path.name.lower() for word in ("protein","peptide","pep","cds"))]
    old_genome=one(root_fastas,"old genomic FASTA at the root of --genome-dir")
    genes_dir=genome_dir.parent/"genes"
    genes=next((genes_dir/name for name in ("genes.GC","genes") if (genes_dir/name).is_file()),None)
    if genes is None:
        raise ValueError(f"Expected genes.GC or genes in sibling directory: {genes_dir}")
    updated_dirs=[path for path in sorted(genome_dir.iterdir()) if path.is_dir()]
    if not updated_dirs: raise ValueError(f"No updated-genome subdirectories in {genome_dir}")
    return old_genome,genes,updated_dirs

def default_output_dir(input_dir):
    """Put remapping results beside other analyses for the standard workspace layout."""
    if input_dir.name == "csv" and input_dir.parent.name == "fitness" and input_dir.parent.parent.name == "analyses":
        return input_dir.parent.parent/"gene_remap"
    if input_dir.name == "csv" and input_dir.parent.name == "fit_analysis":
        return input_dir.parent.parent/"gene_remap_results"
    return input_dir/"gene_remap_results"

def write_results_guide(output, mapping_dirs):
    """Document the remapping outputs and their source evidence."""
    lines=[
        "# Gene remapping results", "",
        "## Start here", "",
        "- `new_locus_info.csv` combines all updated genomes into one review table.",
        "- Its columns give the original locus, description, fitness direction, each updated genome's locus ID and full FAA header, and one PaperBLAST link.",
        "- Each `*_mapping/new_locus_info.csv` file is the detailed mapping table for one updated genome.",
        "- Review rows whose `status` is not `matched_gene`. A `duplicate_protein_assignments.tsv` file is included only when duplicate assignments were found.",
        "", "## Output files and source evidence", "",
        "| Output | Main source | Purpose |",
        "|---|---|---|",
        "| `new_locus_info.csv` | All per-genome mapping tables | Combined review table with one locus-ID and FAA-header pair per updated genome. |",
        "| `*_mapping/new_locus_info.csv` | `all_old_gene_sequences.fasta`, `blast_hits.tsv`, updated-genome GFF3 and protein FASTA | Main locus mapping table for one updated genome, including BLAST quality, mapped gene/protein IDs, annotations, and mapping status. |",
        "| `blast_hits.tsv` | BLAST of `all_old_gene_sequences.fasta` against the updated genomic FASTA | Raw alignment evidence used to choose each best passing hit. |",
        "| `new_locus_sequences.fasta` | Updated genomic FASTA and GFF3 | DNA sequence of each mapped updated gene, or its BLAST-hit interval when no gene overlaps. |",
        "| `new_proteins.faa` | Updated protein FASTA and GFF3 | Protein sequences for mapped genes when a protein could be identified. |",
        "| `positive_gene_sequences.fasta`, `negative_gene_sequences.fasta` | Original genes table and old genome FASTA | Original-gene sequences extracted from the positive and negative input lists. |",
        "| `all_old_gene_sequences.fasta` | The positive and negative extracted FASTA files | Combined BLAST query set. |",
        "",
        "No plots are generated by this pipeline because the results are sequence and locus mappings rather than quantitative measurements.",
    ]
    if mapping_dirs:
        lines.extend(["", "## Updated genomes processed", ""])
        lines.extend(f"- `{path.name}/`" for path in mapping_dirs)
    (output/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(f"Wrote results guide: {output/'README.md'}")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir",type=Path,required=True,help="fitness-analysis csv directory containing top_positive.csv and top_negative.csv")
    p.add_argument("--genome-dir","--genome_dir",dest="genome_dir",type=Path,required=True,
                   help="directory with the old genome FASTA at its root and one updated genome per subdirectory")
    p.add_argument("--genes",type=Path,
                   help="optional genes.GC override; default: sibling genes/genes.GC")
    p.add_argument("--output-dir",type=Path,
                   help="optional output directory; default: workspace gene_remap directory")
    a=p.parse_args(); input_dir=a.input_dir.resolve(); genome_dir=a.genome_dir.resolve()
    output=(a.output_dir or default_output_dir(input_dir)).resolve()
    try:
        if not genome_dir.is_dir(): raise ValueError(f"Genome directory not found: {genome_dir}")
        positive,negative=discover_old_inputs(input_dir)
        old_genome,default_genes,updated_dirs=discover_genomes(genome_dir)
        genes=(a.genes or default_genes).resolve()
        if not genes.is_file(): raise ValueError(f"Genes annotation file not found: {genes}")
        output.mkdir(parents=True,exist_ok=True)
        locs=old.locations(genes); old_sequences=old.fasta(old_genome)
        positive_ids=old.locus_ids(positive); negative_ids=old.locus_ids(negative)
        old.write("positive",positive_ids,locs,old_sequences,"one-based-inclusive",output)
        old.write("negative",negative_ids,locs,old_sequences,"one-based-inclusive",output)
        combined=output/"all_old_gene_sequences.fasta"
        with combined.open("w") as target:
            for path in (output/"positive_gene_sequences.fasta",output/"negative_gene_sequences.fasta"):
                target.write(path.read_text())
        descriptions,categories=gene_annotations(genes); mapping_dirs=[]; mappings=[]
        for directory in updated_dirs:
            genome,gff3,proteins=discover_new_inputs(directory)
            mapping_output=output/f"{directory.name}_mapping"; mapping_output.mkdir(parents=True,exist_ok=True)
            print(f"Old genome: {old_genome.name}\nUpdated genome: {genome.name}\nGFF3: {gff3.name}\nProteins: {proteins.name}")
            result=run_mapping(combined,genome,gff3,proteins,mapping_output,a)
            if result != 0: return result
            enrich_locus_table(mapping_output/"new_locus_info.csv",positive_ids,negative_ids,
                               descriptions,categories,proteins)
            remove_if_empty(mapping_output/"duplicate_protein_assignments.tsv")
            mapping_dirs.append(mapping_output)
            headers,headers_by_base=protein_header_lookup(proteins)
            mappings.append({"directory":mapping_output,"prefix":column_prefix(directory.name),
                             "protein_headers":headers,"protein_headers_by_base":headers_by_base})
        remove_if_empty(output/"positive_missing_locusIds.txt")
        remove_if_empty(output/"negative_missing_locusIds.txt")
        write_summary_table(output,positive_ids,negative_ids,descriptions,mappings)
        write_results_guide(output,mapping_dirs)
        return 0
    except (OSError,ValueError) as error:
        print(f"Error: {error}",file=sys.stderr); return 2

if __name__=="__main__": raise SystemExit(main())
