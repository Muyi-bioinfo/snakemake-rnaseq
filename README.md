# RNA-seq Pipeline

[English](README.md) | [中文](README_CN.md)

[![Snakemake](https://img.shields.io/badge/Snakemake-9+-brightgreen)](https://snakemake.github.io/)
[![Docker](https://img.shields.io/badge/Docker-Supported-blue)](https://www.docker.com/)
[![Apptainer](https://img.shields.io/badge/Apptainer-Supported-orange)](https://apptainer.org/)

A Snakemake-based RNA-seq analysis pipeline covering raw FASTQ QC through differential expression, with full containerization support.

**Toolchain**: FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler (+ gene symbol annotation, GO/KEGG enrichment)

## Features

- **PE/SE auto-adapt** — one `mode` switch in config, all rules adapt
- **Modular rules** — 9 rules in 7 `.smk` modules, each with isolated Conda env
- **Numbered outputs** — `runs/{batch}/results/01_fastqc_raw/` → … → `runs/{batch}/results/07_clusterprofiler/`, execution order reflected in `ls`
- **Three runtime modes** — Conda (dev), Docker (single-machine), Apptainer (HPC)
- **Gene symbol annotation** — DESeq2 output auto-annotated with gene names from the same GTF
- **GO/KEGG enrichment** — clusterProfiler for functional enrichment analysis on DEGs
- **Input validation** — pre-flight check catches missing FASTQ files before execution
- **Sample file** — CSV/TSV/TXT sample table, no hardcoded sample lists in YAML
- **Reproducible** — pinned tool versions across Conda envs and container images
- **Extensible** — add samples or contrasts without touching the Snakefile

📖 **[Snakemake Learning Guide](LEARNING_SNAKEMAKE.md)** — Learn Snakemake core concepts from this project's actual code

## Quick Start

```bash
# 1. Enter project
cd 01.RNA_seq

# 2. Place input files
cp /path/to/sample_*_R1.fastq.gz workflow/data/raw/
cp /path/to/sample_*_R2.fastq.gz workflow/data/raw/
cp genome.fa genes.gtf adapters.fa workflow/data/ref/

# 3. Edit config/samples.tsv — sample names and groups
#    Edit config/config.yaml — contrasts, parameters

# 4. Dry run
snakemake -s workflow/Snakefile -np

# 5. Run
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8
```

## Requirements

| Dependency | Notes |
|---|---|
| [Snakemake](https://snakemake.github.io/) ≥9.0 | Workflow engine |
| [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) | For `--software-deployment-method conda` mode |
| [Docker](https://www.docker.com/) ≥20.04 | For Docker mode |
| [Apptainer](https://apptainer.org/) / [Singularity](https://sylabs.io/) ≥1.0 | For HPC mode |

> Only one of Conda / Docker / Apptainer is required.

## Project Structure

```
.
├── config/
│   ├── config.yaml              # Tool parameters, contrasts, paths
│   └── samples.tsv              # Sample table (TSV/CSV/TXT auto-detect)
├── workflow/                    # Snakemake core
│   ├── Snakefile                #   Orchestrator
│   ├── rules/                   #   7 .smk modules (9 rules)
│   ├── scripts/
│   │   ├── common.py            #   Python helpers (PE/SE-aware I/O)
│   │   ├── common.R             #   R shared utility (parse_arg)
│   │   ├── deseq2.R             #   DESeq2 + diagnostics
│   │   ├── gene2symbol.R        #   Gene ID → Symbol annotation
│   │   └── clusterprofiler.R    #   GO + KEGG enrichment
│   ├── envs/                    #   7 Conda environment files
│   └── data/                    #   Input data
│       ├── raw/                 #     Raw FASTQ
│       └── ref/                 #     Reference genome / GTF / adapters
├── containers/                  # Docker + Apptainer build
│   ├── Dockerfile
│   ├── apptainer.def
│   └── build.sh
├── profile/                     # Snakemake execution profiles
│   ├── docker/
│   ├── apptainer/
│   └── slurm-apptainer/
├── runs/                         # Output (batch-isolated, numbered by DAG layer)
│   └── {batch}/                  #   One subdirectory per run (e.g. 20260701_WT_KO)
│       ├── results/
│       │   ├── 01_fastqc_raw/
│       │   ├── 02_trimmomatic/
│       │   ├── 03_fastqc_trimmed/
│       │   ├── 03_hisat2/
│       │   ├── 04_featurecounts/
│       │   ├── 05_multiqc/
│       │   ├── 06_deseq2/
│       │   └── 07_clusterprofiler/
│       ├── logs/                 #   Per-rule log files
│       └── benchmarks/           #   Per-rule resource usage
├── README.md
└── README_CN.md
```

## Pipeline DAG

```
Layer 0: hisat2_index       One-time genome index
Layer 1: fastqc_raw         Per-sample raw QC
Layer 2: trimmomatic         Per-sample trim  (PE/SE auto)
Layer 3: fastqc_trimmed     Per-sample post-trim QC
Layer 3: hisat2_align       Per-sample alignment (PE/SE auto)
Layer 4: featurecounts      Aggregate gene counts
Layer 5: multiqc            Aggregate QC report
Layer 6: deseq2              Differential expression + gene symbols
Layer 7: clusterprofiler    GO + KEGG enrichment
```

Rules with `{sample}` wildcard run in parallel automatically.

## Configuration

[`config/config.yaml`](config/config.yaml) is the single point of configuration.

### Sample setup

`config/samples.tsv` (TSV format — tab-separated; CSV and TXT also supported by extension):

```tsv
sample	group
WT_1	WT
WT_2	WT
WT_3	WT
KO_1	KO
KO_2	KO
KO_3	KO
```

### PE → SE switch

```yaml
read_pattern:
  mode: "single"             # ← change from "paired"
  single_suffix: ".fastq.gz"
```

No Snakefile or rule changes needed.

### Contrasts

```yaml
deseq2:
  contrasts:
    - name: "KO_vs_WT"
      case: "KO"
      control: "WT"
  padj_threshold: 0.05
  log2fc_threshold: 1.0
```

### ClusterProfiler (enrichment)

```yaml
clusterprofiler:
  org_db: "org.Hs.eg.db"       # species OrgDb (org.Mm.eg.db for mouse)
  kegg_organism: "hsa"         # KEGG code (mmu, rno, ...)
  from_type: "ENSEMBL"         # gene ID type in DEG table → ENTREZID
  pvalue_cutoff: 0.05
  qvalue_cutoff: 0.2
  show_category: 15            # top N terms in dotplots
```

### Switching batches

Change `paths.batch` in `config.yaml` to isolate outputs — all results, logs, and benchmarks go under `runs/{batch}/`:

```yaml
paths:
  batch: "20260715_new_batch"         # ← switch batch
```

All downstream paths auto-follow.

## Usage

### Conda Mode

```bash
snakemake -s workflow/Snakefile -np                    # preview DAG
snakemake -s workflow/Snakefile --dag | dot -Tpng > dag.png
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8       # run
snakemake -s workflow/Snakefile --software-deployment-method conda -j 4 --forcerun trimmomatic
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8 --rerun-incomplete
```

### Docker Mode

```bash
bash containers/build.sh docker                         # build image
bash containers/build.sh test                           # verify tools
snakemake -s workflow/Snakefile --profile profile/docker
```

### Apptainer Mode (HPC)

```bash
bash containers/build.sh all                            # Docker → SIF
snakemake -s workflow/Snakefile --profile profile/apptainer
```

### Slurm + Apptainer

```bash
# Edit profile/slurm-apptainer/config.yaml first (partition, account)
snakemake -s workflow/Snakefile --profile profile/slurm-apptainer
```

## Output Files

### Per-step outputs

| Directory | Key Files |
|---|---|
| `01_fastqc_raw/` | `{sample}_R{1,2}_fastqc.{html,zip}` |
| `02_trimmomatic/` | `{sample}_R{1,2}.trimmed.fastq.gz` |
| `03_fastqc_trimmed/` | `{sample}_R{1,2}.trimmed_fastqc.{html,zip}` |
| `03_hisat2/` | `{sample}.sorted.bam` (+ `.bam.bai`) |
| `04_featurecounts/` | `featurecounts.txt`, `featurecounts.summary.txt` |
| `05_multiqc/` | `multiqc_report.html` |

### DESeq2 output (`06_deseq2/`)

| File | Description |
|---|---|
| `{contrast}_all_results.csv` | All genes with statistics + gene_name column |
| `{contrast}_significant.csv` | Significant DEGs + gene_name column |
| `{contrast}_MA_plot.pdf` | MA plot |
| `{contrast}_volcano_plot.pdf` | Volcano plot (labelled with gene symbols when GTF provided) |
| `PCA_plot.pdf` / `.png` | Sample PCA |
| `sample_distance_heatmap.pdf` | Sample-to-sample distance |
| `DEG_heatmap.pdf` | Top N DEG expression heatmap |

### ClusterProfiler output (`07_clusterprofiler/`)

| File | Description |
|---|---|
| `{contrast}_GO_enrichment.csv` | GO enrichment (BP/CC/MF) results |
| `{contrast}_GO_dotplot.pdf` | GO dotplot — top enriched terms |
| `{contrast}_KEGG_enrichment.csv` | KEGG pathway enrichment results |
| `{contrast}_KEGG_dotplot.pdf` | KEGG pathway dotplot |

## Gene Symbol Annotation

The pipeline automatically annotates DEG output with gene symbols using the same GTF file used for alignment and counting — ensuring complete consistency across the analysis chain.

```bash
# Also usable standalone:
Rscript workflow/scripts/gene2symbol.R \
    --input DEG.csv --gtf genes.gtf --output DEG_anno.csv
```

Implementation: `rtracklayer::import()` (not regex), exact same GTF as `featureCounts` and `hisat2`.

## Container Support

```bash
bash containers/build.sh docker      # Local Docker image
bash containers/build.sh all         # Docker → SIF (recommended for HPC)
bash containers/build.sh pull REGISTRY=docker.io/org  # Remote → SIF
bash containers/build.sh test        # Verify all 11 tools
bash containers/build.sh clean       # Remove images
```

Unified image contains: FastQC 0.12.1, Trimmomatic 0.39, HISAT2 2.2.1, Samtools 1.18, featureCounts 2.0.6, MultiQC 1.21, R 4.3.2 + DESeq2 1.42.0 + rtracklayer 1.62.0 + clusterProfiler 4.10.0 + enrichplot 1.22.0.

## Adding Samples or Contrasts

**New sample** — only `config/samples.tsv`:

```tsv
sample	group
...
NewSample_1	New              # ← append row
```

**New contrast** — only `config/config.yaml`:

```yaml
deseq2:
  contrasts:
    - {name: "New_vs_Ctrl", case: "New", control: "Ctrl"}
```

## Snakemake Features

This pipeline serves as a reference for:

- **Rule definitions** — 9 rules, uniform conda + container + log + benchmark pattern
- **Wildcards** — `{sample}` for 4 per-sample rules, automatic parallelization
- **expand()** — Dynamic target list generation via `get_all_pipeline_targets()`
- **config.yaml** — Centralized parameter management
- **Conda** — Per-rule isolated environments (7 envs)
- **Container** — `docker://` URI shared across all rules
- **Log & Benchmark** — Per-rule tracking
- **DAG** — Explicit input/output contracts, auto-dependency resolution
- **Python functions** — 14 helpers in `common.py`, `os.path.join` output patterns, `*` unpacking
- **PE/SE batching** — Single config switch, all rules adapt automatically
- **Input validation** — Pre-flight FASTQ existence check
- **Sample file** — CSV/TSV/TXT with auto-detected delimiter

## FAQ

**How to add new samples?** — Append rows to `config/samples.tsv`.

**How to change read naming?** — Edit `read_pattern` in `config/config.yaml`.

**How to use single-end reads?** — Set `read_pattern.mode: "single"`. No Snakefile changes.

**How to view benchmarks?** — `cat runs/20260701_WT_KO/benchmarks/trimmomatic/WT_1.txt`

**How to clean restart?** — `rm -rf runs/20260701_WT_KO/ && snakemake -s workflow/Snakefile --software-deployment-method conda -j 8`

**My cluster uses SGE not Slurm?** — Replace `profile/slurm-apptainer/` with SGE config. Snakemake supports SGE, LSF, PBS/Torque, and generic cluster profiles.

## License

For educational and research purposes.
