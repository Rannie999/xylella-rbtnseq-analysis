# Analysis workspace

Create a separate workspace beside this repository, or populate this directory after cloning. Large inputs and generated results are ignored by Git.

```text
workspace/
├── reference/
│   ├── genes/genes_raw             # original gene annotation for MapTnSeq reference
│   ├── genes/genes.GC              # RegionGC.pl output
│   ├── mapping/ML3.tab             # MapTnSeq output for the mutant library
│   ├── genomes/
│   │   ├── TemeculaL.fna            # old genome FASTA matching genes_raw
│   │   ├── TemeculaL_complete/      # updated genome: FASTA, GFF, protein FASTA
│   │   ├── Temecula1/               # updated genome: FASTA, GFF, protein FASTA
│   │   └── 9a5c/                    # updated genome: FASTA, GFF, protein FASTA
│   └── pool/ML3.pool               # DesignRandomPool.pl output
├── barseq_counts/
│   └── EXPERIMENT/*.codes          # BarSeq barcode-count files
├── poolcounts/
│   └── EXPERIMENT/                 # combineBarSeq.pl output
└── runs/
    └── RUN_NAME/
        ├── input/
        │   ├── exps.tsv            # FEBA experiment metadata
        │   ├── sample_name_map.csv         # readable names, aligned to FEBA columns
        │   └── groups.json          # condition-to-replicate mapping
        ├── feba/                    # BarSeqR.pl and RunFEBA.R outputs
        ├── analyses/
        │   ├── barcode_qc/
        │   ├── fitness/
        │   ├── pairwise_ttest/
        │   └── gene_remap/
        └── summaries/
```

`EXPERIMENT` groups `.codes` files sequenced together. `RUN_NAME` identifies the biological experiment and its metadata. Multiple runs can use the same library pool and reference files.
