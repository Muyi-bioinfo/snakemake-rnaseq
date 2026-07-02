# RNA-seq 分析流程

[English](README.md) | [中文](README_CN.md)

[![Snakemake](https://img.shields.io/badge/Snakemake-9+-brightgreen)](https://snakemake.github.io/)
[![Docker](https://img.shields.io/badge/Docker-支持-blue)](https://www.docker.com/)
[![Apptainer](https://img.shields.io/badge/Apptainer-支持-orange)](https://apptainer.org/)

基于 Snakemake 构建的 RNA-seq 上游分析流程，覆盖原始 FASTQ 质控至差异表达基因列表的完整分析链，支持容器化部署。

**分析工具链**：FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler（含基因 Symbol 注释 + GO/KEGG 富集分析）

## 特性

- **PE/SE 自动适配** — config 切换一个 `mode`，全部 rule 自动跟随
- **模块化设计** — 9 个 rule 封装为 7 个 `.smk` 模块，每个有独立 Conda 环境
- **编号输出目录** — `runs/{batch}/results/01_fastqc_raw/` → … → `runs/{batch}/results/07_clusterprofiler/`，`ls` 即见执行顺序
- **三种运行模式** — Conda（开发）、Docker（单机）、Apptainer（HPC 集群）
- **基因 Symbol 注释** — DESeq2 输出自动从同一 GTF 追加 gene symbol 列
- **GO/KEGG 功能富集** — clusterProfiler 对差异基因进行 GO + KEGG 通路富集分析
- **输入校验** — 启动前自动检查 FASTQ 文件是否存在，避免跑到一半报错
- **样本表独立** — 支持 CSV/TSV/TXT 格式样本表，不再硬编码在 YAML 中
- **可复现** — Conda 环境和容器镜像均锁定工具版本
- **易于扩展** — 新增样本或对比仅修改外部文件，零 Snakefile 改动

## 快速开始

```bash
# 1. 进入项目
cd 01.RNA_seq

# 2. 放入测序数据
cp /path/to/sample_*_R1.fastq.gz workflow/data/raw/
cp /path/to/sample_*_R2.fastq.gz workflow/data/raw/

# 3. 放入参考文件
cp genome.fa genes.gtf adapters.fa workflow/data/ref/

# 4. 编辑 config/samples.tsv — 样本名和分组
#    编辑 config/config.yaml — 对比设计、工具参数

# 5. 预览执行计划
snakemake -s workflow/Snakefile -np

# 6. 正式运行
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8
```

## 环境依赖

| 依赖 | 说明 |
|---|---|
| [Snakemake](https://snakemake.github.io/) ≥9.0 | 工作流引擎 |
| [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) | `--software-deployment-method conda` 模式 |
| [Docker](https://www.docker.com/) ≥20.04 | Docker 容器模式 |
| [Apptainer](https://apptainer.org/) / [Singularity](https://sylabs.io/) ≥1.0 | HPC 容器模式 |

> 三者满足其一即可。开发推荐 Conda，生产推荐 Docker/Apptainer。

## 快速导航

| 你想… | 章节 |
|---|---|
| 了解项目结构 | [项目结构](#项目结构) |
| 理解分析流程 | [流程 DAG](#流程-dag) |
| 修改参数/新增样本 | [配置说明](#配置说明) |
| 选择运行方式 | [运行方式](#运行方式) |
| 查看输出结果 | [输出文件说明](#输出文件说明) |
| 在 HPC 上运行 | [Apptainer 模式](#apptainer-模式推荐-hpc) |
| 使用基因 Symbol 注释 | [基因 Symbol 注释](#基因-symbol-注释) |

## 项目结构

```
.
├── config/
│   ├── config.yaml              # 工具参数、对比设计、输出路径
│   └── samples.tsv              # 样本表（TSV/CSV/TXT 自动识别）
├── workflow/                    # Snakemake 工作流核心
│   ├── Snakefile                #   编排器（7 条 include）
│   ├── rules/                   #   7 个 .smk 模块（9 个 rule）
│   │   ├── hisat2.smk           #     hisat2_index + hisat2_align
│   │   ├── fastqc.smk           #     fastqc_raw + fastqc_trimmed
│   │   ├── trimmomatic.smk      #     trimmomatic（PE/SE 自动适配）
│   │   ├── featurecounts.smk    #     featureCounts
│   │   ├── multiqc.smk          #     MultiQC
│   │   ├── deseq2.smk           #     DESeq2（含 gene symbol 注释）
│   │   └── clusterprofiler.smk  #     clusterProfiler（GO + KEGG 富集）
│   ├── scripts/
│   │   ├── common.py            #   Python 辅助函数（PE/SE 感知 I/O）
│   │   ├── common.R             #   R 公共函数（parse_arg）
│   │   ├── deseq2.R             #   DESeq2 差异分析 + 诊断图
│   │   ├── gene2symbol.R        #   基因 ID → Symbol 注释库（rtracklayer）
│   │   └── clusterprofiler.R    #   GO/KEGG 富集分析脚本
│   ├── envs/                    #   7 个 Conda 环境定义
│   └── data/                    #   输入数据
│       ├── raw/                 #     原始 FASTQ
│       └── ref/                 #     参考基因组 / GTF / 接头
├── containers/                  # Docker + Apptainer 构建
│   ├── Dockerfile               #   统合容器镜像
│   ├── apptainer.def          #   Apptainer 定义
│   └── build.sh                 #   双引擎构建脚本
├── profile/                     # Snakemake 执行配置
│   ├── docker/
│   ├── apptainer/
│   └── slurm-apptainer/
├── runs/                         # 输出（按批次隔离，按 DAG 层级编号）
│   └── {batch}/                  #   每次运行一个子目录（如 20260701_WT_KO）
│       ├── results/
│       │   ├── 01_fastqc_raw/    #     FastQC 原始数据
│       │   ├── 02_trimmomatic/   #     Trimmomatic 剪切产物
│       │   ├── 03_fastqc_trimmed/ #    FastQC 剪切后
│       │   ├── 03_hisat2/        #     HISAT2 比对 BAM
│       │   ├── 04_featurecounts/ #     基因计数矩阵
│       │   ├── 05_multiqc/       #     汇总 QC 报告
│       │   ├── 06_deseq2/        #     差异表达结果
│       │   └── 07_clusterprofiler/ #   功能富集分析结果
│       ├── logs/                 #   每个 Rule 的独立日志
│       └── benchmarks/           #   每个 Rule 的资源记录
├── README.md
└── README_CN.md
```

## 流程 DAG

```
Layer 0: hisat2_index      一次性建基因组索引，无 sample wildcard
Layer 1: fastqc_raw        每个样本独立执行，并行
Layer 2: trimmomatic        每个样本独立执行，PE/SE 自动适配
Layer 3: fastqc_trimmed    每个样本独立执行
Layer 3: hisat2_align      每个样本独立执行，PE/SE 自动适配
Layer 4: featurecounts     汇总所有 BAM，单次运行
Layer 5: multiqc           汇总 QC 报告，单次运行
Layer 6: deseq2            读取计数矩阵，单次运行（含 gene symbol 注释）
Layer 7: clusterprofiler   GO + KEGG 功能富集，单次运行，依赖 Layer 6
```

持有 `{sample}` wildcard 的 4 个 rule 由 Snakemake 自动并行。汇总型 rule 等待所有上游样本完成后触发。

## 配置说明

[`config/config.yaml`](config/config.yaml) 是所有可调参数的唯一入口。

### 样本设置

`config/samples.tsv`（TSV 格式，Tab 分隔；也支持 CSV 和 TXT，扩展名自动识别）：

```tsv
sample	group
WT_1	WT
WT_2	WT
WT_3	WT
KO_1	KO
KO_2	KO
KO_3	KO
```

### 双端 → 单端切换

```yaml
read_pattern:
  mode: "single"             # ← 从 "paired" 改为 "single"
  single_suffix: ".fastq.gz"
```

无需修改 Snakefile 或任何 rule 文件。

### 差异比较设计

```yaml
deseq2:
  contrasts:
    - name: "KO_vs_WT"
      case: "KO"
      control: "WT"
  padj_threshold: 0.05
  log2fc_threshold: 1.0
  top_n_genes: 50
```

### clusterProfiler 富集分析

```yaml
clusterprofiler:
  org_db: "org.Hs.eg.db"       # 物种注释库（小鼠用 org.Mm.eg.db）
  kegg_organism: "hsa"         # KEGG 物种代码（mmu, rno, ...）
  from_type: "ENSEMBL"         # DEG 表中的基因 ID 类型 → ENTREZID
  pvalue_cutoff: 0.05
  qvalue_cutoff: 0.2
  show_category: 15            # 气泡图展示的 top N 条目
```

### 切换批次 / 自定义输出目录名

每次运行改 `paths.batch` 即可隔离产出；修改 `paths.outputs.*` 可自定义子目录名，全部代码自动跟随：

```yaml
paths:
  batch: "20260715_new_batch"             # ← 切换批次
  outputs:
    fastqc_raw: "runs/{batch}/results/01_qc_raw"   # ← 改名即可
```

### 新增样本

仅编辑 `config/samples.tsv`：

```tsv
sample	group
...
NewSample_1	New              # ← 追加一行
```

### 新增对比

```yaml
deseq2:
  contrasts:
    - {name: "New_vs_Ctrl", case: "New", control: "Ctrl"}  # ← 添加
```

## 运行方式

### Conda 模式（推荐开发/调试）

```bash
# 预览执行计划
snakemake -s workflow/Snakefile -np

# 生成 DAG 可视化
snakemake -s workflow/Snakefile --dag | dot -Tpng > dag.png

# 正式运行（8 核并行）
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8

# 只运行特定目标
snakemake -s workflow/Snakefile --software-deployment-method conda -j 4 runs/20260701_WT_KO/results/04_featurecounts/featurecounts.txt

# 强制重跑某 rule + 下游
snakemake -s workflow/Snakefile --software-deployment-method conda -j 4 --forcerun trimmomatic

# 中断恢复
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8 --rerun-incomplete
```

> `-s workflow/Snakefile` 是必需的，因为 Snakefile 在 `workflow/` 子目录中（遵循 Snakemake 官方标准布局）。

### Docker 模式

```bash
bash containers/build.sh docker                      # 构建镜像
bash containers/build.sh test                        # 验证工具
snakemake -s workflow/Snakefile --profile profile/docker
```

### Apptainer 模式（推荐 HPC）

```bash
bash containers/build.sh all                         # 构建 Docker → 转 SIF
snakemake -s workflow/Snakefile --profile profile/apptainer
```

### Slurm + Apptainer（多节点 HPC）

使用前需修改 `profile/slurm-apptainer/config.yaml` 的分区和账户：

```bash
snakemake -s workflow/Snakefile --profile profile/slurm-apptainer
```

## 输出文件说明

### 各步骤输出

| 目录 | 关键文件 |
|---|---|
| `01_fastqc_raw/` | `{sample}_R{1,2}_fastqc.{html,zip}` |
| `02_trimmomatic/` | `{sample}_R{1,2}.trimmed.fastq.gz` |
| `03_fastqc_trimmed/` | `{sample}_R{1,2}.trimmed_fastqc.{html,zip}` |
| `03_hisat2/` | `{sample}.sorted.bam` + `.bam.bai` |
| `04_featurecounts/` | `featurecounts.txt`、`featurecounts.summary.txt` |
| `05_multiqc/` | `multiqc_report.html` |

### DESeq2 输出（`06_deseq2/`）

| 文件 | 内容 |
|---|---|
| `{contrast}_all_results.csv` | 全部基因统计表，含 gene_name 列 |
| `{contrast}_significant.csv` | 显著差异基因列表，含 gene_name 列 |
| `{contrast}_MA_plot.pdf` | MA 图（平均表达量 vs. log2FC） |
| `{contrast}_volcano_plot.pdf` | 火山图（有 GTF 时标注 gene symbol） |
| `PCA_plot.pdf` / `.png` | 样本 PCA 降维图 |
| `sample_distance_heatmap.pdf` | 样本间表达距离矩阵 |
| `DEG_heatmap.pdf` | Top N 差异基因表达热图 |

### clusterProfiler 输出（`07_clusterprofiler/`）

| 文件 | 内容 |
|---|---|
| `{contrast}_GO_enrichment.csv` | GO 富集分析结果（BP/CC/MF） |
| `{contrast}_GO_dotplot.pdf` | GO 富集气泡图 |
| `{contrast}_KEGG_enrichment.csv` | KEGG 通路富集结果 |
| `{contrast}_KEGG_dotplot.pdf` | KEGG 通路气泡图 |

## 基因 Symbol 注释

流程自动使用与比对和计数同一份 GTF 文件，为 DEG 输出追加基因 Symbol 列，确保上游到下游注释版本一致。

```bash
# 也可独立使用
Rscript workflow/scripts/gene2symbol.R \
    --input DEG.csv --gtf genes.gtf --output DEG_anno.csv
```

实现方式：`rtracklayer::import()` 标准解析 GTF（非正则），与 `featureCounts` 和 `hisat2` 用同一文件。

## 容器化支持

```bash
bash containers/build.sh docker      # 本地 Docker 镜像
bash containers/build.sh all         # Docker → SIF（HPC 推荐）
bash containers/build.sh pull REGISTRY=docker.io/org  # 远程拉取 → SIF
bash containers/build.sh test        # 验证镜像内全部 11 个工具
bash containers/build.sh clean       # 清理
```

统合镜像包含：FastQC 0.12.1、Trimmomatic 0.39、HISAT2 2.2.1、Samtools 1.18、featureCounts 2.0.6、MultiQC 1.21、R 4.3.2 + DESeq2 1.42.0 + rtracklayer 1.62.0 + clusterProfiler 4.10.0 + enrichplot 1.22.0。

## Snakemake 特性展示

| 特性 | 应用 | 说明 |
|---|---|---|
| **Rule 定义** | 9 个 rule | conda + container + log + benchmark 统一模式 |
| **Wildcards** | `{sample}` | 4 个 per-sample rule 自动并行 |
| **expand()** | `get_all_pipeline_targets()` | 顶层目标列表动态生成 |
| **config.yaml** | 全局配置 | 参数/对比/输出路径/容器 集中管理 |
| **Conda** | 7 个 env yaml | 每个 rule 的隔离环境 |
| **Container** | `docker://` URI | 全部 rule 共用统合镜像 |
| **日志 & Benchmark** | `runs/{batch}/logs/` `runs/{batch}/benchmarks/` | 每个 rule 独立追踪，按批次隔离 |
| **DAG** | input→output 链 | 显式依赖，自动推导并行 |
| **Python 函数** | `common.py` 14 函数 | I/O、模式、目标生成、校验、样本表加载 |
| **PE/SE 批处理** | config 一个开关 | 全部受影响 rule 自动适配 |
| **输入校验** | 启动前检查 | 缺失 FASTQ 提前报错 |
| **样本表** | CSV/TSV/TXT | 扩展名自动识别分隔符 |

## 常用维护操作

### 查看 Benchmark

```bash
grep -E "s|cpu_time|max_rss" runs/20260701_WT_KO/benchmarks/trimmomatic/WT_1.txt
```

### 清理并重启

```bash
rm -rf runs/20260701_WT_KO/
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8
```

### 日志查看

```bash
cat runs/20260701_WT_KO/logs/trimmomatic/WT_1.log        # 特定样本
tail -f runs/20260701_WT_KO/logs/hisat2_align/*.log      # 实时查看
```

## 工具版本清单

| 工具 | 版本 | Conda 环境 |
|---|---|---|
| FastQC | 0.12.1 | `envs/fastqc.yaml` |
| Trimmomatic | 0.39 | `envs/trimmomatic.yaml` |
| HISAT2 | 2.2.1 | `envs/hisat2.yaml` |
| Samtools | 1.18 | `envs/hisat2.yaml` |
| featureCounts | 2.0.6 | `envs/featurecounts.yaml` |
| MultiQC | 1.21 | `envs/multiqc.yaml` |
| R | 4.3.2 | `envs/deseq2.yaml` |
| DESeq2 | 1.42.0 | `envs/deseq2.yaml` |
| rtracklayer | 1.62.0 | `envs/deseq2.yaml` |
| clusterProfiler | 4.10.0 | `envs/clusterprofiler.yaml` |
| enrichplot | 1.22.0 | `envs/clusterprofiler.yaml` |

## 常见问题

**如何新增样本？** — 在 `config/samples.tsv` 中追加一行即可。

**如何改读段命名格式？** — 编辑 `config/config.yaml` 中的 `read_pattern` 段。

**如何用单端数据？** — 设 `read_pattern.mode: "single"`，无需改 Snakefile。

**如何查看 Benchmark？** — `cat runs/20260701_WT_KO/benchmarks/trimmomatic/WT_1.txt`

**如何彻底重跑？** — `rm -rf runs/20260701_WT_KO/ && snakemake -s workflow/Snakefile --software-deployment-method conda -j 8`

**集群用 SGE 不是 Slurm？** — 替换 `profile/slurm-apptainer/` 为 SGE 配置。Snakemake 支持 SGE、LSF、PBS/Torque 及通用集群配置。

## 许可证

本项目仅供教育和研究目的使用。
