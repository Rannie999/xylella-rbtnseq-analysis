# Detailed workflow

This document records the directory layout and command order used for the Xylella RB-TnSeq analysis. It is intended for reproducing a run from a mapped mutant library and BarSeq `.codes` files. The commands use a concrete example workspace, count batch, and run name; replace them for another experiment.

## Separate code from data

Keep this repository as the code directory and create a data workspace beside it. The [workspace template](workspace/README.md) defines the expected layout. This separation permits version control of commands and scripts without committing raw reads, genome files, count matrices, or results.

```bash
git clone <repository-url> xylella-rbtnseq-analysis
cd xylella-rbtnseq-analysis
```

Run every command below from the workspace root, `workspace/`. The examples process the count batch `exp4_B` and create the run `exp4_B_092126`.

## Requirements

Create and activate a Bioconda environment with Python 3.10, the Python packages in `requirements.txt`, Perl with DBI, R, and BLAST+:

```bash
conda create -n rbtnseq -c conda-forge -c bioconda \
  python=3.10 matplotlib numpy pandas scikit-learn scipy seaborn statsmodels \
  perl perl-dbi r-base blast
conda activate rbtnseq
```

FEBA commands require the Perl modules used by the bundled FEBA code. Gene remapping requires `makeblastdb` and `blastn` on `PATH`.

## Define the R path

Activate the environment that contains R, then use `command -v Rscript` to find its executable. Define `RSCRIPT` once in the workspace shell before running the pool or FEBA steps:

```bash
export RSCRIPT="$(command -v Rscript)"
echo "$RSCRIPT"
"$RSCRIPT" --version
```

Directory structure is shown below. Reference files are shared across runs; each run directory contains its own metadata, FEBA outputs, analyses, and summaries.

```text
xylella-rbtnseq-workspace/
├── reference/
│   ├── genes/
│   │   ├── genes                 # original annotation matching the ML3 library
│   │   └── genes.GC                  # RegionGC.pl output
│   ├── mapping/
│   │   └── Xylella_ML3.tab           # MapTnSeq output
│   ├── genomes/
│   │   ├── TemeculaL.fna               # old genome FASTA matching genes_raw
│   │   ├── TemeculaL_complete/         # updated genome: FASTA, GFF, protein FASTA
│   │   ├── Temecula1/                  # other genome: FASTA, GFF, protein FASTA
│   │   └── 9a5c/                       # other genome: FASTA, GFF, protein FASTA
│   └── pool/
│       └── ML3.pool                  # DesignRandomPool.pl output
├── barseq_counts/
│   └── EXPERIMENT/
│       └── *.codes                   # BarSeq count files from one count batch
├── poolcounts/
│   └── EXPERIMENT/
│       ├── EXPERIMENT.poolcount      # combined barcode-count matrix
│       ├── EXPERIMENT.colsum
│       └── combineBarSeq.log
└── runs/
    └── RUN_NAME/
        ├── input/
        │   ├── exps.tsv              # FEBA sample and Time0 assignment table
        │   ├── sample_name_map.csv           # FEBA-column to readable-name mapping
        │   └── groups.json            # condition-to-replicate mapping
        ├── feba/
        │   ├── all.poolcount
        │   ├── fit_logratios.tab
        │   ├── fit_t.tab
        │   ├── fit_quality.tab
        │   └── *.log
        ├── analyses/
        │   ├── barcode_qc/
        │   ├── fitness/
        │   ├── pairwise_ttest/
        │   └── gene_remap/
        └── summaries/
            ├── analysis_summary.md
            ├── barcode_gene_summary.csv
            └── fitness_summary.csv
```

## Inputs and derived reference files

### Library reference inputs

These source files are shared across all experiments made from the same ML3 mutant pool:

| File | Location in workspace | Source | Purpose |
| --- | --- | --- | --- |
| `Xylella_ML3.tab` | `reference/mapping/` | MapTnSeq output | Barcode-to-insertion mapping evidence. |
| `genes` | `reference/genes/` | Matching genome annotation | Gene coordinates for pool and FEBA processing. |
| `TemeculaL.fna` | `reference/genomes/` | Genome assembly used for mapping | Required to calculate gene GC content. |

### Generated library reference files

Generate these once from the matching source inputs, then reuse them across runs from the same library:

| Generated file | Location in workspace | Created from | Purpose |
| --- | --- | --- | --- |
| `ML3.pool` | `reference/pool/` | `Xylella_ML3.tab` and `genes` with `DesignRandomPool.pl` | One accepted insertion assignment per barcode. |
| `genes.GC` | `reference/genes/` | `genes` and `TemeculaL.fna` with `RegionGC.pl` | Annotation with `GC` and `nTA`, required by `BarSeqR.pl`. |

All annotation, pool, and FASTA files must use the same scaffold identifiers.

### Run-specific input

| File | Location in workspace | Created from | Purpose |
| --- | --- | --- | --- |
| `*.codes` | `barseq_counts/EXPERIMENT/` | BarSeq barcode-count output | Per-sample barcode counts that `combineBarSeq.pl` combines into one poolcount matrix. |
| `exps.tsv` | `runs/RUN_NAME/input/` | Curated experiment metadata | Maps each poolcount column to a sample and assigns its Time0 baseline. |
| `sample_name_map.csv` | `runs/RUN_NAME/input/` | Curated sample-name mapping | Maps FEBA sample-column names to readable names for Barcode QC and fitness analysis. |
| `groups.json` | `runs/RUN_NAME/input/` | Curated condition metadata | Defines the experimental condition and replicate samples for fitness and pairwise analyses. |

These files are not shared across experiments: make an `exps.tsv`, `sample_name_map.csv`, and `groups.json` file for each `RUN_NAME`, and retain them with that run's results.

## Prepare inputs for runs

### 1. Generate the library pool

`DesignRandomPool.pl` converts MapTnSeq output into one reliable barcode-to-insertion row per mutant. Reuse the resulting pool for every BarSeq experiment made from that ML3 library.

**Inputs:** `reference/mapping/Xylella_ML3.tab` and `reference/genes/genes`.

**Output:** `reference/pool/ML3.pool`, the reusable barcode-to-insertion pool definition, and `logs/DesignRandomPool.log`

```bash
/usr/bin/perl /path/to/feba/bin/DesignRandomPool.pl \
  -pool /path/to/pool_file.pool \
  -genes /path/to/genes \
  /path/to/mapping_output.tab \
  > /path/to/DesignRandomPool.log 2>&1
```
example:

```bash
/usr/bin/perl ../feba/bin/DesignRandomPool.pl \
  -pool reference/pool/ML3.pool \
  -genes reference/genes/genes \
  reference/mapping/Xylella_ML3.tab \
  > logs/DesignRandomPool.log 2>&1
```

`DesignRandomPool.pl` invokes `PoolStats.R` after writing the pool. If `PoolStats.R` cannot find `Rscript`, activate the R environment; the pool itself may still be complete while diagnostic `.hit`, `.unhit`, and `.surprise` files are absent.

### 2. Generate the GC annotation once using the matching genome assembly

`BarSeqR.pl` needs a gene table with a `GC` column. Generate it from the original annotation and the matching genome FASTA (`TemeculaL.fna`):

**Inputs:** `reference/genes/TemeculaL.fna` and `reference/genes/genes`.

**Output:** `reference/genes/genes.GC`, the gene annotation augmented with GC content and `nTA`.

```bash
/usr/bin/perl /path/to/feba/bin/RegionGC.pl \
  /path/to/genome.fna \
  /path/to/genes \
  > /path/to/genes.GC
```
example:
```bash
/usr/bin/perl ../feba/bin/RegionGC.pl \
  reference/genes/TemeculaL.fna \
  reference/genes/genes \
  > reference/genes/genes.GC
```
> **Note:** `scaffoldId` values in `TemeculaL.fna`, `genes` and `EXPERIMENT.poolcount` (see step 1 in Step-by-Step Python Workflow below) must be consistent.

### 3. Prepare `exps.tsv`, `sample_name_map.csv`, and `groups.json`

`BarSeqR.pl` uses `exps.tsv` to map each count-table column to a sample and choose its Time0 baseline. The essential columns are:

```text
SetName    Index    Description    Date_pool_expt_started    Group
```

- `Index` must exactly match a sample column in `EXPERIMENT.poolcount`.
- Set `Group` to `Time0` for baseline samples.
- Set the **same** values for all samples intended to share one baseline. `SetName` in `exps.tsv` must match the prefix in `EXPERIMENT.poolcount`.FEBA first assigns Time0 samples that share both `SetName` and `Date_pool_expt_started` with the experiment. If no exact match exists, it first falls back to Time0 samples from another lane on the same date, then to a Time0 sample from the same set on another date. 

`sample_name_map.csv` supplies `old_colnames` and `new_colnames`: `old_colnames` must match the FEBA count or fitness column name, and `new_colnames` is the readable sample name used in downstream outputs. `groups.json` assigns each experimental condition to its sample names. Its names must use the same names that appear after applying `sample_name_map.csv`; Time0 samples are excluded from fitness-condition groups.

## Run FEBA pipeline

### Quick-run shell scripts

The `workflow/` directory provides shell-script wrappers that derive standard locations from their directory arguments. Run them from the workspace directory. They create expected output directories, stop when a required input is missing, and print result and log paths after a successful run. They provide a shorter way to run the workflow; the detailed commands are in the `Step-by-Step Python Workflow`.

```bash
# Combine .codes files for one sequencing experiment (step 1).
bash path/to/workflow/run_combineBarSeq.sh path/to/codes_directory

# Prepare FEBA inputs and all.poolcount (step 2).
# The run input directory must contain exps.tsv.
bash path/to/workflow/run_BarSeqR.sh path/to/poolcounts/exp4_B path/to/runs/exp4_B_092126/input

# Estimate fitness with FEBA (step 3). Set RSCRIPT if Rscript is not on PATH.
bash path/to/workflow/run_runFEBA.sh path/to/experiment_run_dir

# Run downstream analyses (step 4).
# The run input directory must contain exps.tsv, sample_name_map.csv, groups.json.
# Run individual analyses: run_barcode_qc.sh, run_fitness.sh, run_pairwise_ttest.sh, run_gene_remap.sh
bash path/to/workflow/SCRIPT path/to/experiment_run_dir

# Or run all four downstream analyses in sequence.
bash path/to/workflow/run_downstream_analyses.sh path/to/experiment_run_dir
```
Example commands:
```bash
bash ../workflow/run_combineBarSeq.sh barseq_counts/exp4_B

bash ../workflow/run_BarSeqR.sh poolcounts/exp4_B runs/exp4_B_092126/input

bash ../workflow/run_runFEBA.sh runs/exp4_B_092126

bash ../workflow/run_barcode_qc.sh runs/exp4_B_092126
bash ../workflow/run_fitness.sh runs/exp4_B_092126
bash ../workflow/run_pairwise_ttest.sh runs/exp4_B_092126
bash ../workflow/run_gene_remap.sh runs/exp4_B_092126

bash ../workflow/run_downstream_analyses.sh runs/exp4_B_092126
```


- `run_downstream_analyses.sh` runs Barcode QC, fitness analysis, pairwise tests, and gene remapping in that order. It accepts an optional second argument for the FEBA t-like threshold (default to 2.5) for `run_fitness.sh`.
- `run_fitness.sh` accepts the same optional threshold, for example `bash ../workflow/run_fitness.sh runs/exp4_B_092126 3`. 
- `run_gene_remap.sh` uses `runs/exp4_B_092126/analyses/fitness/csv/` and `reference/genomes/` automatically.

## Step-by-Step Python Workflow

### 1. Combine BarSeq counts

The raw data `.codes` files contain per-sample barcode counts. `combineBarSeq.pl` reverse-complements each barcode and matches it to the library pool. It produces one matrix including all samples in one experiment.

**Inputs:** `barseq_counts/exp4_B/*.codes` and `reference/pool/ML3.pool`.

**Outputs:** `poolcounts/exp4_B/exp4_B.poolcount`, `poolcounts/exp4_B/exp4_B.colsum`, and `poolcounts/exp4_B/combineBarSeq.log`.

`combineBarSeq.pl` usage:
```bash
/usr/bin/perl /path/to/feba/bin/combineBarSeq.pl \
  -all \
  /path/to/output/output_prefix \
  /path/to/pool_file.pool \
  /path/to/*.codes \
  > /path/to/combineBarSeq.log 2>&1
```
example:

```bash
mkdir -p poolcounts/exp4_B
/usr/bin/perl ../feba/bin/combineBarSeq.pl \
  -all \
  poolcounts/exp4_B/exp4_B \
  reference/pool/ML3.pool \
  barseq_counts/exp4_B/*.codes \
  > poolcounts/exp4_B/combineBarSeq.log 2>&1
```

This creates `EXPERIMENT.poolcount`, `EXPERIMENT.colsum`, and `combineBarSeq.log`. The wrapper uses `-all`, which retains samples that fail FEBA's default count or pool-match filters for later QC. It does not create counts for barcodes absent from the supplied pool.

### 2. Build `all.poolcount`

Run `BarSeqR.pl` once for each analysis run:

This step links the combined barcode counts to the pool's insertion and gene assignments, then prepares the FEBA input directory.

**Inputs:** `poolcounts/exp4_B/exp4_B.poolcount`, `runs/exp4_B_092126/input/exps.tsv`, `reference/genes/genes.GC`, and `reference/pool/ML3.pool`.

**Outputs:** `runs/exp4_B_092126/feba/all.poolcount`, `runs/exp4_B_092126/feba/BarSeqR.log`, `runs/exp4_B_092126/feba/pool`, and `runs/exp4_B_092126/feba/genes`. The `pool` and `genes` files are direct copies of the supplied `ML3.pool` and `genes.GC` files.

`BarSeqR.pl` usage:
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
example:
```bash
mkdir -p runs/exp4_B_092126/feba
/usr/bin/perl ../feba/bin/BarSeqR.pl \
  -org Xylella \
  -indir poolcounts/exp4_B/ \
  -exps runs/exp4_B_092126/input/exps.tsv \
  -genesfile reference/genes/genes.GC \
  -poolfile reference/pool/ML3.pool \
  -outdir runs/exp4_B_092126/feba \
  -noR \
  > runs/exp4_B_092126/feba/BarSeqR.log 2>&1
```

The wrapper creates `runs/RUN_NAME/feba/all.poolcount` and writes `BarSeqR.log`. It uses `-noR` so the R fitness step remains explicit. FEBA's optional `-metadir` gives media and compound labels; it does not change barcode counts. Before a production run, add `-test` to the `BarSeqR.pl` command if you want FEBA to validate the count-table and metadata relationship without producing final results.

### 3. Estimate fitness with FEBA

`RunFEBA.R` calculates strain and gene fitness by comparing each experimental sample with the Time0 baseline assigned through `exps.tsv`.

**Inputs:** `runs/exp4_B_092126/feba/all.poolcount`, the generated FEBA input files in that directory, and the bundled FEBA code.

**Outputs:** `fit_logratios.tab`, `fit_t.tab`, `strain_fit.tab`, and `fit_quality.tab`, etc. under `runs/exp4_B_092126/feba/`.

`RunFEBA.R` usage:
```bash
/path/to/Rscript \
  /path/to/feba/bin/RunFEBA.R \
  orgname \
  /path/to/feba_output_directory \
  /path/to/feba_code \
  > /path/to/runFEBA.log 2>&1
```
example:
```bash
Rscript \
  ../feba/bin/RunFEBA.R \
  Xylella \
  runs/exp4_B_092126/feba \
  ../feba \
  > runs/exp4_B_092126/feba/runFEBA.log 2>&1
```

The bundled `feba/lib/FEBA.R` is customized for this Xylella dataset. It lowers the minimum genes per scaffold from 10 to 3, lowers the minimum `cor12` quality threshold from 0.10 to 0.05, and issues a warning rather than stopping when fewer than 100 `genesUsed12` genes are available. It also writes `d1_gN` and `d2_gN` diagnostics, prints quality metrics, and uses serial `lapply()` instead of `mclapply()`. Review the quality outputs before interpreting results because these changes retain lower-quality evidence than upstream defaults.

### 4. Run downstream analyses

Downstream analyses were performed using custom scripts outside the FEBA repository, including:

* **Barcode QC**, which evaluates barcode counts and insertion coverage.
* **Fitness analysis**, which summarizes gene-fitness effects across replicates.
* **Pairwise testing**, which compares gene-fitness effects between biological conditions.
* **Gene mapping**, which maps genes to the other Xylella genomes.

#### 4.1. Barcode QC

**Inputs:** `runs/exp4_B_092126/feba/all.poolcount` and the run-specific `runs/exp4_B_092126/input/sample_name_map.csv`.

`sample_name_map.csv` is an optional comma-separated name-mapping table. `barcode_qc.py` requires these exact columns when it is supplied:

| Column | Purpose |
| --- | --- |
| `old_colnames` | Exact sample column name in `all.poolcount`. |
| `new_colnames` | Readable replacement name used in QC tables and plots. |

Additional columns are ignored by this script. Omit `--name-map` to retain the original `all.poolcount` sample names.

**Outputs:** `runs/exp4_B_092126/analyses/barcode_qc/`.
The output includes sample-level read depth, detected barcodes and genes, low-count fractions, Time0 retention, insertion coverage, and plots. CSV files retain complete values even when plot labels are abbreviated.

`barcode_qc.py` usage:
```bash
python scripts/barcode_qc.py \
  --input /path/to/all.poolcount \
  --name-map /path/to/sample_name_map.csv \
  --output /path/to/barcode_qc
```
example:
```bash
python ../scripts/barcode_qc.py \
  --input runs/exp4_B_092126/feba/all.poolcount \
  --name-map runs/exp4_B_092126/input/sample_name_map.csv \
  --output runs/exp4_B_092126/analyses/barcode_qc
```

#### 4.2. Fitness analysis

**Inputs:** `runs/exp4_B_092126/feba/fit_logratios.tab`, `runs/exp4_B_092126/feba/fit_t.tab`, and the run-specific `runs/exp4_B_092126/input/sample_name_map.csv` and `runs/exp4_B_092126/input/groups.json` files.

`groups.json` defines the replicate samples for each biological condition. Each key is the condition label used in results, and its value is a list of readable sample names that exactly match the `new_colnames` values in `sample_name_map.csv`:
```json
{
  "condition_A": ["condition_A_rep1", "condition_A_rep2"],
  "condition_B": ["condition_B_rep1", "condition_B_rep2"]
}
```

**Outputs:** `runs/exp4_B_092126/analyses/fitness/`.

`fit_analysis.py` usage:
```bash
python scripts/fit_analysis.py \
  --input /path/to/feba_output_directory \
  --output /path/to/fitness \
  --name-map /path/to/sample_name_map.csv \
  --groups /path/to/groups.json \
  --threshold THRESHOLD
```
example:
```bash
python ../scripts/fit_analysis.py \
  --input runs/exp4_B_092126/feba \
  --output runs/exp4_B_092126/analyses/fitness \
  --name-map runs/exp4_B_092126/input/sample_name_map.csv \
  --groups runs/exp4_B_092126/input/groups.json \
  --threshold 2.5
```

#### 4.3. Pairwise condition tests

The formal analysis runs Welch tests for every valid gene in each condition pair. Benjamini-Hochberg q-values are calculated separately per condition comparison across all valid gene tests. The `--filtered` option runs the test using results filtered based on t-like values.

**Inputs:** `runs/exp4_B_092126/feba/fit_logratios.tab`, `runs/exp4_B_092126/feba/fit_t.tab`, `runs/exp4_B_092126/feba/analyses/fitness/csv/fit_filtered.csv`, and the run-specific `runs/exp4_B_092126/input/sample_name_map.csv` and `runs/exp4_B_092126/input/groups.json` files.

**Outputs:** `runs/exp4_B_092126/analyses/pairwise_ttest/`.

`pairwise_ttest.py` usage:
```bash
python path/to/scripts/pairwise_ttest.py \
  --input /path/to/feba_output_directory \
  --name-map /path/to/sample_name_map.csv \
  --groups /path/to/groups.json \
  --filtered /path/to/fitness/csv/fit_filtered.csv \
  --output /path/to/pairwise_ttest
```
example:
```bash
python ../scripts/pairwise_ttest.py \
  --input runs/exp4_B_092126/feba \
  --name-map runs/exp4_B_092126/input/sample_name_map.csv \
  --groups runs/exp4_B_092126/input/groups.json \
  --filtered runs/exp4_B_092126/analyses/fitness/csv/fit_filtered.csv \
  --output runs/exp4_B_092126/analyses/pairwise_ttest
```

#### 4.4. Gene remapping

This optional step extracts candidate sequences from the fitness-analysis output and maps them to the other Xylella genome annotations.

**Inputs:** `runs/exp4_B_092126/analyses/fitness/csv/top_positive.csv` or `runs/exp4_B_092126/analyses/fitness/csv/top_negative.csv`, `reference/genes/genes.GC`, the original genome FASTA, and the updated-genome directories, containing one genomic FASTA, one GFF/GFF3, and one protein FASTA.

**Outputs:** `runs/exp4_B_092126/analyses/gene_remap/`, including `new_locus_info.csv` and sequence-mapping evidence.

`gene_remap.py` usage:
```bash
python scripts/gene_remap.py \
  --input-dir /path/to/fitness_csv_directory \
  --genome-dir /path/to/genomes \
  --genes /path/to/genes.GC \
  --output-dir /path/to/gene_remap_results
```
example:
```bash
python ../scripts/gene_remap.py \
  --input-dir runs/exp4_B_092126/analyses/fitness/csv/ \
  --genome-dir reference/genomes \
  --genes reference/genes/genes.GC \
  --output-dir runs/exp4_B_092126/analyses/gene_remap
```

`--genome-dir` contains the old genome FASTA directly in that directory and one subdirectory for each other genome.

The root `new_locus_info.csv` is the combined review table. Each `*_mapping/` subdirectory contains the detailed table and BLAST evidence for one updated genome. Its generated README explains the mapping tables, FASTA files, mapping status, and duplicate-protein report. This step generates mapping and sequence evidence rather than quantitative plots.

#### 4.5. Summarize a completed FEBA run

**Inputs:** `runs/exp4_B_092126/` for one run, or  `runs/` when summarizing all runs.

**Outputs:** `runs/exp4_B_092126/summaries/` for one run, or `runs/multi_run_summary/` when summarizing all runs.

`summarize_feba_run.py` usage:
```bash
python scripts/summarize_feba_run.py --run_dir /path/to/run_directory
```
example:
```bash
# Summarize exp4_B
python ../scripts/summarize_feba_run.py --run_dir runs/exp4_B_092126

# Summarize every completed run under runs/.
python ../scripts/summarize_feba_run.py --run_dir runs/
```

This writes `analysis_summary.md`, `barcode_gene_summary.csv`, and `fitness_summary.csv`. A single-run summary is written to that run's `summaries/` directory. A `runs/` summary is written to `runs/multi_run_summary/`. Both CSV tables have one row per run, with the run-directory name in the `Experiment` column. The barcode/gene table reports barcode, gene, read-depth, insertion, intergenic, and Time0-retention metrics. The fitness table reports represented genes, possible and non-missing values, large-effect and t-like passing values, quality-passing samples, and mean absolute fitness and t-like values.


## Output locations and review order

| Step | Main output | First review |
| --- | --- | --- |
| Pool generation | `reference/pool/ML3.pool` | Mapping log and optional PoolStats summaries. |
| Count combining | `poolcounts/EXPERIMENT/EXPERIMENT.poolcount` | `EXPERIMENT.colsum` and log. |
| BarSeqR | `runs/RUN_NAME/feba/all.poolcount` | Barcode and genic insertion counts in `BarSeqR.log`. |
| RunFEBA | `runs/RUN_NAME/feba/fit_quality.tab` | `fit_quality.tab`, correlation plots, and `genesUsed12`. |
| Barcode QC | `analyses/barcode_qc/` | Sample read depth, barcode detection, and Time0 retention. |
| Fitness analysis | `analyses/fitness/` | Its README, `condition_summary.csv`, and high-confidence calls. |
| Pairwise tests | `analyses/pairwise_ttest/` | Its README and `comparison_summary.csv`. |
| Gene remapping | `analyses/gene_remap/` | `new_locus_info.csv` and its README. |

For a reproducible record, retain each run's `input/` files and logs, record the Git commit with `git rev-parse HEAD`, and do not overwrite an existing `RUN_NAME` after modifying inputs or thresholds.
