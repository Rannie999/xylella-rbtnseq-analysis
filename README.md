# Xylella RB-TnSeq analysis

Reproducible quality control and downstream analysis for Xylella RB-TnSeq and BarSeq experiments. The bundled FEBA code performs barcode mapping and gene-fitness estimation; the custom scripts analyze FEBA outputs, compare conditions, and remap candidate genes to other Xylella genomes.

The `feba/` directory contains the FEBA code used by this workflow. It originates from the [Berkeley Lab FEBA repository](https://bitbucket.org/berkeleylab/feba/src/master/); see [feba/LICENSE](feba/LICENSE) and [feba/CITATION](feba/CITATION). The bundled `feba/lib/FEBA.R` includes dataset-specific changes described below.

For a complete reproducible workspace layout, metadata rules, Bash wrappers, command order, and output review guide, read [WORKFLOW.md](WORKFLOW.md).

## Analysis flow

```text
MapTnSeq output -> DesignRandomPool.pl -> library pool
                                            |
BarSeq .codes -> combineBarSeq.pl -----------+
                                            |
                                     combined poolcount
                                            |
genes.GC + exps.tsv -> BarSeqR.pl -> all.poolcount
                                            |
                                        RunFEBA.R
                                            |
                fit_logratios.tab + fit_t.tab + fit_quality.tab
                                            |
                  QC, fitness analysis, pairwise tests, gene remapping
```

## Requirements

Create and activate a Bioconda environment with Python 3.10, the Python packages in `requirements.txt`, Perl with DBI, R, and BLAST+:

```bash
conda create -n rbtnseq -c conda-forge -c bioconda \
  python=3.10 matplotlib numpy pandas scikit-learn scipy seaborn statsmodels \
  perl perl-dbi r-base blast
conda activate rbtnseq
```

FEBA commands require the Perl modules used by the bundled FEBA code. Gene remapping requires `makeblastdb` and `blastn` on `PATH`.

## Prepare inputs for runs

### 1. Build the ML3 pool definition

`DesignRandomPool.pl` converts MapTnSeq output into one reliable barcode-to-insertion row per mutant. Reuse the resulting pool for every BarSeq experiment made from that ML3 library.

With the FEBA defaults, a barcode needs at least 10 good mapping reads, 75% support for its preferred location, and an eight-fold margin over the next-best location. The pool, gene table, and reference genome must use identical scaffold identifiers.

```bash
/usr/bin/perl /path/to/feba/bin/DesignRandomPool.pl \
  -pool /path/to/pool_file.pool \
  -genes /path/to/genes \
  /path/to/mapping_output.tab \
  > /path/to/DesignRandomPool.log 2>&1
```


`DesignRandomPool.pl` writes the pool before it invokes `PoolStats.R`. If `PoolStats.R` cannot execute, the pool can still be valid, but the `.hit`, `.unhit`, and `.surprise` summaries will be missing. Activate the environment that provides `Rscript` before running this step.

### 2. Add GC content to the gene annotation

`BarSeqR.pl` needs a gene table with a `GC` column. `RegionGC.pl` reads the genome assembly that was used to map the mutant library and the original gene-coordinate table, calculates GC content (adds `GC` column) and the number of possible insertion sites (add `nTA` column) for each gene, and writes the augmented table to standard output. Redirect that output to `genes.GC`.

```bash
/usr/bin/perl /path/to/feba/bin/RegionGC.pl \
  /path/to/genome.fna \
  /path/to/genes \
  > /path/to/genes.GC
```

### 3. Prepare `exps.tsv`, `sample_name_map.csv`, and `groups.json`

`exps.tsv` maps poolcount columns to samples and assigns Time0 baselines. Its essential columns are `SetName`, `Index`, `Description`, `Date_pool_expt_started`, and `Group`. `Index` must match the poolcount sample column, `Group = Time0` identifies baseline samples, and `SetName` must match the poolcount filename prefix. FEBA first uses Time0 samples with matching `SetName` and `Date_pool_expt_started`; if needed, it falls back to another lane on the same date, then the same set on another date.

`sample_name_map.csv` contains `old_colnames` and `new_colnames`, mapping FEBA column names to readable downstream-analysis names. `groups.json` maps each biological condition to the readable replicate names in `new_colnames`; omit Time0 samples from these condition groups.

These three files have separate roles: FEBA uses `exps.tsv` to locate counts and choose Time0 controls, while the custom Python analyses use `sample_name_map.csv` and `groups.json` to rename samples and identify biological replicates.

## Run FEBA pipeline

### 1. Combine BarSeq barcode counts

`combineBarSeq.pl` matches reverse-complement barcodes from `.codes` files to the pool. One run produces one shared matrix with barcodes as rows and samples as columns.

The first positional argument is the output prefix. FEBA appends `.poolcount` and `.colsum` to it. The second positional argument is the reusable pool; all remaining arguments are one or more `.codes` files from the same sequencing experiment.

```bash
/usr/bin/perl /path/to/feba/bin/combineBarSeq.pl \
  -all \
  /path/to/output_prefix \
  /path/to/pool_file.pool \
  /path/to/*.codes \
  > /path/to/combineBarSeq.log 2>&1
```

This writes `output_prefix.poolcount` and `output_prefix.colsum`. `-all` retains samples that do not meet FEBA's default read-depth or pool-matching filters so that they can be inspected during QC. It does not add barcodes absent from the pool.

### 2. Create `all.poolcount`

`BarSeqR.pl` creates the complete barcode-by-sample input table for FEBA. `EXPERIMENT.poolcount` is the observed-count matrix from the current BarSeq samples. It includes insertion-coordinate columns copied by `combineBarSeq.pl`, but it is not the authoritative mutant list: it can omit accepted pool barcodes that have zero counts in every current sample.

`-poolfile` supplies that complete reference list of accepted mutants. `BarSeqR.pl` uses it to retain every pool barcode, match barcode coordinates consistently, and copy the exact reference to the FEBA output directory. It then uses `genes.GC` to assign each insertion to a `locusId`, adds all sample-count columns from the poolcount file, and writes `all.poolcount`. The resulting file therefore includes every accepted barcode, its insertion location and gene assignment, and its observed count in each sample.

`-indir` contains `EXPERIMENT.poolcount`, `-exps` supplies the run-specific sample and Time0 metadata, and `-outdir` is the run's `feba/` directory. `-noR` deliberately stops after preparing `pool`, `genes`, `exps`, and `all.poolcount` so the R fitness calculation can be run and logged separately.

```bash
/usr/bin/perl /path/to/feba/bin/BarSeqR.pl \
  -org Xylella \
  -indir /path/to/poolcount_directory \
  -exps /path/to/exps.tsv \
  -genesfile /path/to/genes.GC \
  -poolfile /path/to/pool_file.pool \
  -outdir /path/to/feba_output_directory \
  -noR \
  > /path/to/BarSeqR.log 2>&1
```

`-metadir` is optional. It supplies FEBA's media and compound labels and does not change barcode counts. Add `-test` before a production run to validate poolcount names, sample columns, and metadata without writing final results.

### 3. Calculate FEBA fitness

#### Dataset-specific `FEBA.R` changes

The bundled `feba/lib/FEBA.R` differs from the upstream FEBA version to retain more data for this Xylella dataset. It relaxes three default quality filters:

| Criterion | Upstream FEBA | Bundled code |
| --- | ---: | ---: |
| Minimum genes per scaffold | 10 | 3 |
| Minimum `cor12` quality threshold | 0.10 | 0.05 |
| Fewer than 100 `genesUsed12` genes | Stop analysis | Issue a warning and continue |

The customized code also writes `d1_gN` and `d2_gN` half-gene count diagnostics, prints pseudovariance and quality metrics to the log, and uses serial `lapply()` in place of parallel `mclapply()`. These diagnostic and execution changes do not alter the quality thresholds.

Relaxing the checks can retain additional experiments and candidate genes, but it also permits lower-quality evidence than the upstream defaults. Review `fit_quality.tab`, `fit_quality_cor12.pdf`, and the reported `genesUsed12` values before interpreting results.

`RunFEBA.R` reads `all.poolcount`, `genes`, `pool`, and `exps` from the FEBA output directory. It assigns the Time0 baseline for each sample, calculates strain and gene fitness, performs FEBA quality checks, and saves result tables and diagnostic plots in that same directory. The three positional arguments are the organism label, the prepared FEBA output directory, and the bundled FEBA code directory.

```bash
/path/to/Rscript \
  /path/to/feba/bin/RunFEBA.R \
  orgname \
  /path/to/feba_output_directory \
  /path/to/feba_code \
  > /path/to/runFEBA.log 2>&1
```

The key outputs are `fit_logratios.tab`, `fit_t.tab`, `strain_fit.tab`, and `fit_quality.tab`. Negative fitness means mutants in that gene declined relative to their matched Time0 baseline; positive fitness means they increased.

## 4. Analyze FEBA outputs

### 4.1. Barcode QC

`barcode_qc.py` evaluates count-matrix (`all.poolcount`) quality before interpreting fitness. It summarizes read depth, pool detection, genic and intergenic insertion coverage, and retention of Time0 barcodes for each sample. It does not change FEBA input or fitness files.

```bash
python scripts/barcode_qc.py \
  --input /path/to/all.poolcount \
  --name-map /path/to/sample_name_map.csv \
  --output /path/to/barcode_qc
```

`sample_name_map.csv` is an optional comma-separated name-mapping table. `barcode_qc.py` requires these exact columns when it is supplied:

| Column | Purpose |
| --- | --- |
| `old_colnames` | Exact sample column name in `all.poolcount`. |
| `new_colnames` | Readable replacement name used in QC tables and plots. |

Additional columns are ignored by this script. Omit `--name-map` to retain the original `all.poolcount` sample names.

The output includes sample-level read depth, detected barcodes and genes, low-count fractions, Time0 retention, insertion coverage, and plots. CSV files retain complete values even when plot labels are abbreviated.

### 4.2. Fitness analysis

`fit_analysis.py` reads FEBA gene fitness and t-like statistics, removes Time0 columns from downstream comparisons, renames experimental columns, and summarizes replicate agreement by condition. `--threshold` is the absolute t-like cutoff used when classifying supporting replicate evidence; it does not recalculate FEBA fitness values.

`--input` is the FEBA output directory and must contain both `fit_logratios.tab` and `fit_t.tab`.

`groups.json` defines the replicate samples for each biological condition. Each key is the condition label used in results, and its value is a list of readable sample names that exactly match the `new_colnames` values in `sample_name_map.csv`:

```json
{
  "condition_A": ["condition_A_rep1", "condition_A_rep2"],
  "condition_B": ["condition_B_rep1", "condition_B_rep2"]
}
```

```bash
python scripts/fit_analysis.py \
  --input /path/to/feba_output_directory \
  --output /path/to/fitness \
  --name-map /path/to/sample_name_map.csv \
  --groups /path/to/groups.json \
  --threshold 2.5
```

Use only experimental samples in `groups.json`; the file excludes names beginning with `T0-` as Time0 baselines. 

The output directory contains a README, summary and candidate CSV files under `csv/`, and QC and fitness plots under `plots/`. `csv/high_confidence_genes.csv` contains the primary replicated findings; the output README identifies the source CSV for each plot.

### 4.3. Pairwise condition tests

`pairwise_ttest.py` compares every valid gene between each pair of conditions defined in `groups.json`. It uses the complete FEBA fitness matrix for formal Welch tests and Benjamini–Hochberg q-values. `--filtered` provides a separate exploratory analysis of candidates selected by `fit_analysis.py`. Their q-values are conditional on preselection and should not be interpreted as formal FDR-controlled results.

`--input` is the FEBA output directory and must contain both `fit_logratios.tab` and `fit_t.tab`. The optional `--filtered` file is `analyses/fitness/csv/fit_filtered.csv` produced by `fit_analysis.py`.

```bash
python scripts/pairwise_ttest.py \
  --input /path/to/feba_output_directory \
  --name-map /path/to/sample_name_map.csv \
  --groups /path/to/groups.json \
  --filtered /path/to/fitness/csv/fit_filtered.csv \
  --output /path/to/pairwise_ttest
```

### 4.4. Gene remapping

Candidate loci from `top_positive.csv` and `top_negative.csv`, which are outputs of `fit_analysis.py`, can be remapped to an updated genome:

`gene_remap.py` extracts the original DNA sequences for those candidate loci, aligns them against each updated genome, and uses its GFF/GFF3 and protein FASTA to identify the corresponding locus and protein. The root result table combines assignments across updated genomes; the per-genome directories retain the detailed BLAST evidence.

`--genome-dir` should contain the old genome FASTA directly in that directory and one subdirectory for each updated genome. Each updated-genome subdirectory must contain one genomic FASTA, one GFF/GFF3 file, and one protein FASTA. 

```bash
python scripts/gene_remap.py \
  --input-dir /path/to/fitness_csv_directory \
  --genome-dir /path/to/genomes \
  --genes /path/to/genes.GC \
  --output-dir /path/to/gene_remap_results
```

The root `new_locus_info.csv` is the combined review table. It contains `query_locusId`, `description`, `fitness_change`, a `locusId` and full `faa_header` for every provided genome, and `paperblast_url`. Each `*_mapping/` directory contains the detailed `new_locus_info.csv` table and BLAST evidence for one provided genome. Its generated README explains the mapping tables, FASTA files, and mapping status. This step generates mapping and sequence evidence rather than quantitative plots.

### 4.5. Summarize a completed FEBA run

`summarize_feba_run.py` reads completed FEBA output and creates compact coverage and fitness-quality summaries. Give it one run directory to summarize that run, or the parent `runs/` directory to discover and compare every immediate subdirectory containing `feba/`.

```bash
python scripts/summarize_feba_run.py --run_dir /path/to/run_directory
```

Pass a `runs/` directory to summarize every completed run it contains. The script writes `analysis_summary.md`, `barcode_gene_summary.csv`, and `fitness_summary.csv`; both CSV files have one row per run, with the run-directory name in the `Experiment` column. They cover barcode/gene representation, read depth, Time0 retention, fitness-effect rates, t-like passing rates, and FEBA quality status.

## Reproducibility

Record the command, threshold, metadata, replicate groups, FEBA version, and Git commit for each analysis:

```bash
git rev-parse HEAD
python scripts/fit_analysis.py --help
python scripts/pairwise_ttest.py --help
```

Raw reads, poolcounts, FEBA outputs, generated analyses, FASTA files, and figures are intentionally ignored by Git. Keep them in the analysis workspace or an archival data repository.
