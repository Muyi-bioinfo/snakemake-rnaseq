# CLAUDE.md

> Snakemake RNA-seq 上游分析流程 — AI 辅助开发上下文

## 项目概述

基于 Snakemake 的 RNA-seq 双端/单端上游分析流程，覆盖原始 FASTQ → 差异基因列表的完整分析链。

**工具链**: `FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler`  
**输入**: FASTQ (PE `{sample}_R1/R2.fastq.gz` 或 SE `{sample}.fastq.gz`)  
**输出**: 差异基因 CSV + PCA/MA/Volcano/Heatmap + GO/KEGG 富集分析 + 基因 Symbol 注释

## 目录结构

```
01.RNA_seq/
├── config/
│   └── config.yaml              # 唯一配置入口 — 样本/参数/对比/输出路径
├── workflow/                    # Snakemake 工作流核心
│   ├── Snakefile                #   编排器 (7 include)
│   ├── rules/                   #   7 个 .smk 模块 (9 rules)
│   │   ├── hisat2.smk           #     hisat2_index + hisat2_align
│   │   ├── fastqc.smk           #     fastqc_raw + fastqc_trimmed
│   │   ├── trimmomatic.smk      #     trimmomatic (PE/SE 自动适配)
│   │   ├── featurecounts.smk    #     featureCounts
│   │   ├── multiqc.smk          #     MultiQC
│   │   ├── deseq2.smk           #     DESeq2 (+ gene symbol)
│   │   └── clusterprofiler.smk  #     clusterProfiler (GO + KEGG)
│   ├── scripts/
│   │   ├── common.py            #   Python 辅助函数
│   │   ├── common.R             #   R 公共函数 (parse_arg)
│   │   ├── deseq2.R             #   DESeq2 差异分析主脚本
│   │   ├── gene2symbol.R        #   基因 ID→Symbol 注释库
│   │   └── clusterprofiler.R    #   GO/KEGG 富集分析脚本
│   ├── envs/                    #   7 个 Conda 环境定义
│   └── data/                    #   输入数据
│       ├── raw/                 #     原始 FASTQ
│       └── ref/                 #     参考基因组 / GTF / 接头
├── containers/                  # 容器化支持
│   ├── Dockerfile               #   统合镜像 (全部工具)
│   ├── apptainer.def          #   Apptainer 定义
│   └── build.sh                 #   双引擎构建脚本
├── profile/                     # Snakemake 执行配置
│   ├── docker/config.yaml
│   ├── apptainer/config.yaml
│   └── slurm-apptainer/config.yaml
├── runs/                         # 运行时产出 (按批次隔离)
│   └── {batch}/                  #   每次运行一个子目录 (如 20260701_WT_KO)
│       ├── results/              #     分析产出 (编号目录)
│       │   ├── 01_fastqc_raw/
│       │   ├── 02_trimmomatic/
│       │   ├── 03_fastqc_trimmed/
│       │   ├── 03_hisat2/
│       │   ├── 04_featurecounts/
│       │   ├── 05_multiqc/
│       │   ├── 06_deseq2/
│       │   └── 07_clusterprofiler/
│       ├── logs/                 #     每个 Rule 的独立日志
│       └── benchmarks/           #     每个 Rule 的运行时记录
├── README.md
├── README_CN.md
└── CLAUDE.md
```

## DAG 拓扑

```
Layer 0: hisat2_index      一次性建索引，无 wildcard
Layer 1: fastqc_raw        每个 {sample} 一次
Layer 2: trimmomatic       每个 {sample} 一次，依赖 Layer 1
Layer 3: fastqc_trimmed    每个 {sample} 一次，依赖 Layer 2
Layer 3: hisat2_align      每个 {sample} 一次，依赖 Layer 2 + Layer 0
Layer 4: featurecounts     汇总所有 BAM，单次运行
Layer 5: multiqc           汇总所有 QC 报告，单次运行
Layer 6: deseq2            读取计数矩阵，单次运行，依赖 Layer 4
Layer 7: clusterprofiler   GO + KEGG 富集，单次运行，依赖 Layer 6
```

Layer 1-3 的四条规则持有 `{sample}` wildcard，Snakemake 自动并行执行所有样本。

## 架构设计要点

### Rule 统一模板

每个 rule 遵循相同结构：`conda:` (开发) + `container:` (生产) + `log:` + `benchmark:`。所有 shell 块首行加 `set -euo pipefail`。输出路径通过三层体系避免硬编码：

```
config/config.yaml  →  paths.outputs.{step} + paths.logs + paths.benchmarks   (模板，含 {batch} 占位符)
workflow/Snakefile  →  加载后 .format(batch=...) 统一替换 → OUT_* / LOG_ROOT / BENCH_ROOT 常量
workflow/rules/*    →  os.path.join(OUT_*/LOG_ROOT/BENCH_ROOT, ...) 引用 (Python)
                     →  {config[paths][outputs][xxx]} 引用 (shell)
```

批次切换只需改 `config.yaml` 中 `paths.batch` 一个值。

### Conda 环境约定

所有 `envs/*.yaml` 中 channel 顺序固定为 `conda-forge → bioconda → defaults`，符合 Bioconda 官方推荐，避免 conda-forge 的新版本被 bioconda 旧版本覆盖。

### PE/SE 模式切换

`config.yaml` 中 `read_pattern.mode: "paired"|"single"` 控制全局。三个受影响的 rule (`fastqc.smk`, `trimmomatic.smk`, `hisat2.smk`) 通过顶层 Python `if/else` 在 rule 定义级切换 PE/SE，input key 和 output 文件数量不同。`rule all` 通过 `get_all_pipeline_targets()` 动态生成目标列表。不需要平行 rule。

### Snakefile → rules 拆分

`Snakefile` 仅保留 `import`、`configfile`、`SAMPLES`、`OUT_*` 常量、`rule all`、7 条 `include:`。每个 `rules/*.smk` 头部有块注释列出 PE/SE 两种模式下的实际 input/output 路径。

### common.py 函数职责

| 函数 | 用途 | 被引用位置 |
|------|------|-----------|
| `get_read_input_list()` | 返回读段文件列表 (PE→2, SE→1) | `fastqc.smk` |
| `get_read_inputs()` | 返回命名读段 dict (PE→{r1,r2}, SE→{read}) | `trimmomatic.smk` |
| `get_trimmed_input_list()` | 剪切后的文件列表 | `fastqc.smk` |
| `get_hisat2_input()` | 比对输入 dict | `hisat2.smk` |
| `get_fastqc_raw_outputs()` | 原始 FastQC 输出列表 | `fastqc.smk` |
| `get_fastqc_trimmed_outputs()` | 剪切后 FastQC 输出列表 | `fastqc.smk` |
| `get_trimmomatic_outputs()` | Trimmomatic 输出列表 | `trimmomatic.smk` |
| `get_all_pipeline_targets()` | rule all 完整目标列表 | `Snakefile` |
| `get_deseq2_outputs()` | DESeq2 输出文件列表 | `Snakefile`, `deseq2.smk` |
| `get_deseq2_groups_json()` | 分组信息 JSON | `deseq2.smk` |
| `get_deseq2_contrasts_json()` | 对比设计 JSON | `deseq2.smk` |
| `get_clusterprofiler_outputs()` | clusterProfiler 输出文件列表 | `Snakefile`, `clusterprofiler.smk` |
| `load_samples_table()` | 从 CSV/TSV 加载样本表 | `get_samples`, `get_deseq2_groups_json` |
| `validate_input_files()` | 启动前校验 FASTQ 存在 | `Snakefile` |

核心路径取值：`_out(config, step)` → `config["paths"]["outputs"][step]`

### R 脚本协作

```
deseq2.smk  --gtf {input.gtf} →
  deseq2.R:
    if (has_gtf) {
      source("workflow/scripts/gene2symbol.R")
      anno <- load_gtf_annotation(gtf)     # rtracklayer::import()
      res_df <- add_gene_symbol(res_df, anno)
      sig    <- add_gene_symbol(sig, anno)
    }

clusterprofiler.smk  --sig_files {params.sig_files_str} →
  clusterprofiler.R:
    bitr(gene_ids, fromType = from_type, toType = "ENTREZID", OrgDb = org_db)
    enrichGO(gene = genes_entrez, OrgDb = org_db, ont = ...)
    enrichKEGG(gene = genes_entrez, organism = kegg_org)
```

`gene2symbol.R` 支持两种调用方式：
- **独立运行**: `Rscript gene2symbol.R --input DEG.csv --gtf genes.gtf`
- **被 source()**: `load_gtf_annotation(gtf)` + `add_gene_symbol(df, anno)`

### R 脚本日志规范

三份 R 脚本统一使用 `[模块名]` 前缀消息格式（`message("[DESeq2] ...")` / `[clusterProfiler]` / `[gene2symbol]`），便于 `grep` 从 Snakemake 日志中过滤特定模块输出。

## 配置约定

[config/config.yaml](config/config.yaml) 是唯一配置入口。关键段：

| 段 | 内容 |
|---|---|
| `paths` | 数据目录、`batch` 批次标识、`outputs.*`/`logs`/`benchmarks` 子目录 (含 `{batch}` 占位符) |
| `container` | 镜像名、registry、SIF 文件路径 |
| `samples_file` | 样本表路径 (CSV/TSV/TXT，扩展名自动识别) |
| `read_pattern` | `mode` (paired/single)、PE 后缀、SE 后缀 |
| `reference` | 基因组 FASTA、GTF、HISAT2 索引前缀 |
| `fastqc` / `trimmomatic` / `hisat2` / `featurecounts` | 各工具参数 |
| `deseq2` | 对比设计、显著性阈值、VST 参数 |
| `clusterprofiler` | 物种 OrgDb、KEGG 代码、富集阈值 |

### 新增样本

仅修改 `config/samples.tsv`:

```tsv
sample	group
WT_1	WT
NewSample_1	New              # ← 添加
```

文件为 TSV 格式（Tab 分隔），也可使用 CSV（逗号分隔）或 TXT——扩展名自动识别分隔符。

### 新增对比

```yaml
deseq2:
  contrasts:
    - name: "New_vs_Ctrl"    # ← 添加
      case: "New"
      control: "Ctrl"
```

### 切换 PE/SE

```yaml
read_pattern:
  mode: "single"             # ← 改这里
```

无需修改 Snakefile 或任何 rule。

### 切换批次 / 自定义输出目录名

每次运行改 `paths.batch` 即可隔离产出；修改 `paths.outputs.*` 可自定义子目录名，所有代码自动跟随：

```yaml
paths:
  batch: "20260715_new_batch"             # ← 切换批次
  outputs:
    fastqc_raw: "runs/{batch}/results/01_qc_raw"   # ← 改名
```

## 运行命令

```bash
# 从项目根目录执行

# preview DAG
snakemake -s workflow/Snakefile -np

# Conda 模式 (开发)
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8

# Docker 模式
bash containers/build.sh docker
snakemake -s workflow/Snakefile --profile profile/docker

# Apptainer 模式 (HPC)
bash containers/build.sh all
snakemake -s workflow/Snakefile --profile profile/apptainer

# Slurm + Apptainer (多节点)
snakemake -s workflow/Snakefile --profile profile/slurm-apptainer

# 单独运行某个目标（路径需含 batch，以实际 config.yaml 中 paths.batch 值为准）
snakemake -s workflow/Snakefile --software-deployment-method conda -j 4 runs/20260701_WT_KO/results/04_featurecounts/featurecounts.txt

# 强制重跑某 rule + 下游
snakemake -s workflow/Snakefile --software-deployment-method conda -j 4 --forcerun trimmomatic

# 中断恢复
snakemake -s workflow/Snakefile --software-deployment-method conda -j 8 --rerun-incomplete
```

> `-s workflow/Snakefile` 是必需的，因为 Snakefile 不在项目根目录（遵循 Snakemake 官方标准布局）。

## 工具版本

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

## 容器化

`containers/build.sh` 是容器操作的统一入口：

```bash
bash containers/build.sh docker         # 构建 Docker 镜像
bash containers/build.sh all            # Docker + 转 SIF (推荐)
bash containers/build.sh test           # 验证镜像内所有工具
bash containers/build.sh clean          # 清理
```

镜像基于 `continuumio/miniconda3:24.1.2-0`，所有工具通过 conda 安装。

## Snakemake 注意事项

1. **Shell 健壮性**: 所有 shell 块首行必须加 `set -euo pipefail`，确保管道错误（如 `hisat2 | samtools sort`）不会静默吞噬失败
2. **Shell 中 JSON 参数**: 传给 R 脚本的 JSON 用单引号包裹 `'{params.groups_json}'`，防止 shell 解析 JSON 双引号
3. **touch() 哨兵文件**: `hisat2_index` 使用 `touch()` 创建哨兵——HISAT2 生成 8 个 `.ht2` 文件，哨兵文件避免逐个检查
4. **GTF 作为 input**: DESeq2 的 GTF 放在 `input:` 而非 `params:`，shell 中引用 `{input.gtf}`，确保 GTF 更新后 Snakemake 自动重跑
5. **`*` 输出展开**: `rule all` 和 `rule deseq2` 使用 `*get_deseq2_outputs(config)` 将 Python 列表展开为多个 input/output 项
6. **output 路径拼接**: 统一使用 `os.path.join(OUT_XXX, "{sample}.html")` 模式，Snakemake 9.x 禁止 output 函数/λ
7. **input 函数 λ 包装**: input 函数需 `lambda wildcards: func(wildcards, config)` 闭包传递 config，因函数在 common.py 模块中无法访问 Snakemake 全局变量
8. **PE/SE rule 级切换**: `fastqc.smk`, `trimmomatic.smk`, `hisat2.smk` 用 `if config[...]` 在 rule 定义层切换 PE/SE，避免 shell 模板中引用不存在的 input key
9. **configfile 路径**: Snakefile 在 `workflow/` 内，configfile 是 `"config/config.yaml"` (Snakemake 9.x 相对 CWD 解析，非 Snakefile 目录)
10. **import 路径**: 用 `Path(_scripts).is_dir()` 存在性检测代替 `__file__`，Snakemake 9.x 中 `__file__` 指向自身模块
11. **include 语义**: Snakemake 的 `include:` 是文本级插入，Snakefile 中定义的 `SAMPLES`、`OUT_*`、`CONTAINER_URI`、Python import 在子模块中直接可用
12. **conda 和 container 互斥**: `--software-deployment-method conda` 时不激活 container；容器模式设置 `software-deployment-method: [apptainer]`（镜像已包含全部工具）
13. **批次路径格式化**: Snakefile 在 `configfile:` 后立即对 `paths.outputs`、`paths.logs`、`paths.benchmarks` 做 `.format(batch=...)` 统一替换 `{batch}` 占位符。`.smk` 规则中 `log:`/`benchmark:` 使用 `os.path.join(LOG_ROOT/BENCH_ROOT, "...")`，shell 中通过 `{config[paths][xxx]}` 引用已格式化的完整路径
14. **mkdir -p 不重复创建 log/bench 子目录**: Snakemake 自动为 `log:` 和 `benchmark:` 创建父目录，shell 中 `mkdir -p` 只需创建 output 目录即可

## 不属于本项目的文件

- `RNA_seq_chatgpt/`、`RNA_seq_deepseek/`、`RNA_seq_gemini/` — 其他参考实现
- `workflow/scripts/__pycache__/` — Python 字节码缓存
