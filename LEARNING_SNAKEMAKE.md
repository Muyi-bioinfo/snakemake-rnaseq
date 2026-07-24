# LEARNING_SNAKEMAKE.md

> 系统化 Snakemake 学习指南 — 从基础概念到企业级部署
>
> 以本项目 RNA-seq 上游分析 pipeline 为实战案例，每章结合真实代码讲解。
>
> 四部分结构：**快速开始 → 核心概念 → 实战进阶 → 附录**

---

## 前置约定

- 终端命令以 `$` 开头，代码块标注语言（```python / ```yaml / ```bash / ```r）
- 示例代码来自项目实际文件，代码块外标注可点击的文件链接（链接在代码块后面，不在前面）
- `{batch}` 是 config 占位符（Snakefile 中 `.format(batch=...)` 替换），`{sample}` 是 Snakemake wildcard（运行时自动匹配）

---

# 第一部分：快速开始

> 目标：跑起来，建立 DAG 思维，能写简单 Rule，会调试

## 01. Snakemake 介绍与 DAG 原理

### 什么是 Workflow

生物信息学分析通常由多个步骤串联而成：质控 → 剪切 → 比对 → 定量 → 差异分析 → 富集分析。每一步依赖前一步的输出。Workflow 引擎将这种**步骤依赖**和**数据流向**形式化描述，自动按正确顺序执行。

### Snakemake 原理

Snakemake 是**基于规则的文件驱动模式**工作流引擎。核心思想：

1. 你定义 **Rule**（规则）：给定 input 文件，产出 output 文件
2. Snakemake 从最终目标**逆向推导**依赖关系，构建 DAG（有向无环图）
3. 自动并行执行无依赖关系的任务

### 为什么选择 Snakemake

| 特性 | Snakemake | Nextflow | WDL | CWL |
|------|-----------|----------|-----|-----|
| 语法基础 | Python 扩展 | Groovy DSL | WDL 方言 | YAML/JSON |
| 学习曲线 | 平缓（Python 用户友好） | 中等 | 中等 | 偏陡 |
| 容器支持 | Docker/Apptainer (原生) | Docker/Apptainer/Sarus | Docker | Docker |
| HPC 支持 | Slurm/SGE/LSF/PBS | 各平台 Executor | Cromwell 后端 | HPC 后端 |
| 社区生态 | 中等 | 大 | 中等 | 中等 |
| 检查点恢复 | 内置（`--rerun-incomplete`） | 内置（resume） | 依赖后端 | 依赖后端 |

Snakemake 的优势在于：Python 原生集成（input/output 可直接写 Python 表达式）、声明式规则定义、自动并行化、内置 Conda 集成。

### DAG 思想：从 output 逆向推导 input

Snakemake 不关心"先跑什么后跑什么"——你只告诉它**最终要什么文件**，它自己算出执行顺序。

```
给定最终目标: results/06_deseq2/KO_vs_WT_significant.csv

Snakemake 逆向推导:
  KO_vs_WT_significant.csv
    ← deseq2 rule (需要 featurecounts.txt + groups + contrasts)
      ← featurecounts rule (需要所有样本的 .sorted.bam + GTF)
        ← hisat2_align rule × N个样本 (需要 trimmed FASTQ + 索引)
          ← trimmomatic rule × N个样本 (需要 raw FASTQ)
            ← fastqc_raw rule × N个样本 (需要 raw FASTQ)  ← 兄弟依赖
          ← hisat2_index rule (需要 genome.fa)             ← 索引是独立分支
```

### DAG 深度与并行度

同一层级无依赖的 job 自动并行。本项目中，Layer 1-3 的规则持有 `{sample}` wildcard，Snakemake 自动为每个样本生成独立的 job，用 `-j` 控制并行数：

```
Layer 0: hisat2_index          (1 job, 无 wildcard)
Layer 1: fastqc_raw            (N jobs, {sample} × N)
Layer 2: trimmomatic           (N jobs, {sample} × N)
Layer 3: fastqc_trimmed        (N jobs, {sample} × N)
Layer 3: hisat2_align          (N jobs, {sample} × N)
Layer 4: featurecounts         (1 job, 汇总所有 BAM)
Layer 5: multiqc               (1 job, 汇总 QC)
Layer 6: deseq2                (1 job)
Layer 7: clusterprofiler       (1 job)
```

### `include:` 文本级插入语义

Snakemake 的 `include:` 是**文本级插入**（类似 C 的 `#include`），不是模块导入。Snakefile 中定义的 `SAMPLES`、`OUT_*`、`CONTAINER_URI`、Python import 在子模块 `.smk` 中直接可用。

```python
# Snakefile:126-147 — 7 条 include 加载规则模块
# Layer 0 + 3b: 基因组索引 + 序列比对
include: "rules/hisat2.smk"

# Layer 1 + 3a: 原始质控 + 剪切后质控
include: "rules/fastqc.smk"

# Layer 2: 去接头 + 质量剪切
include: "rules/trimmomatic.smk"

# Layer 4: 基因表达定量
include: "rules/featurecounts.smk"

# Layer 5: 汇总质控报告
include: "rules/multiqc.smk"

# Layer 6: 差异表达分析
include: "rules/deseq2.smk"

# Layer 7: 基因功能富集分析
include: "rules/clusterprofiler.smk"
```
> 🔗 [Snakefile:126-147](workflow/Snakefile#L126-L147) — 7 条 `include:` 按 DAG 层级排列，顺序不影响执行

### 可视化命令

```bash
$ snakemake -s workflow/Snakefile --dag | dot -Tpdf > dag.pdf
```

> 🔗 [README_CN.md#流程-dag](README_CN.md#流程-dag) — 流程 DAG 图解
> 🔗 [Snakefile:1-34](workflow/Snakefile#L1-L34) — DAG 概览 + 模块映射注释
> 🔗 [Snakefile:114-123](workflow/Snakefile#L114-L123) — `rule all` 是 DAG 的根节点

### ⭐ 最佳实践

- **从 `rule all` 开始设计**：先想清楚最终要什么文件，再反向写 rule
- **每个 rule 只做一件事**：一个 rule = 一个工具的一个运行模式
- **include 顺序不影响执行**：Snakemake 按 DAG 依赖推导，不按文本顺序
- **把 DAG 画出来**：`--dag | dot -Tpdf` 是理解和调试的最佳工具

### ⚠ 常见错误

- **在 Snakefile 中写死路径**：应通过 config + Python 函数动态生成
- **把 include 当作 Python import**：include 是文本插入，不是模块导入
- **忘记 rule all 必须声明所有最终目标**：漏掉一个文件 = 对应的 rule 永远不会执行

### 💡 面试参考 → 附录 D-1, D-2

---

## 02. Rule 编写

### Rule 结构

```python
rule rule_name:
    input:   ...     # 输入文件（依赖声明）
    output:  ...     # 输出文件（产出声明）
    params:  ...     # 非文件参数
    conda:   ...     # Conda 环境
    container: ...   # 容器镜像
    log:     ...     # 日志文件
    benchmark: ...   # 性能基准
    threads: ...     # 线程数
    resources: ...   # 资源限制
    shell:   "..."   # Shell 命令（三选一）
    # run:   "..."   # Python 代码（三选一）
    # script: "..."  # 外部脚本（三选一）
```

### `input:` — 文件依赖声明

支持多种形式：

```python
# 简单文件
input: "data/raw/sample1.fastq.gz"

# 命名 key（在 shell 中通过 {input.key} 引用）
input:
    fasta = "data/ref/genome.fa",

# 列表
input:
    bams = ["bam1.bam", "bam2.bam"],

# Python 函数（动态生成，接收 wildcards 参数）
input:
    lambda wildcards: get_read_input_list(wildcards, config),
```
> 🔗 [fastqc.smk:13](workflow/rules/fastqc.smk#L13) — `lambda wildcards:` input 函数

### `output:` — 文件产出声明

**这是 Snakemake 最关键的指令**。Snakemake 根据 output 文件是否存在、时间戳是否新于 input，决定是否需要运行该 rule。

```python
output:
    bam = os.path.join(OUT_HISAT2, "{sample}.sorted.bam"),
    bai = os.path.join(OUT_HISAT2, "{sample}.sorted.bam.bai"),
```
> 🔗 [hisat2.smk:39-41](workflow/rules/hisat2.smk#L39-L41) — 命名 output key，`{sample}` wildcard 自动匹配

### `shell:` vs `run:` vs `script:`

| 指令 | 用途 | 典型场景 |
|------|------|---------|
| `shell:` | 执行 shell 命令 | 调用命令行工具（FastQC, HISAT2, Trimmomatic） |
| `run:` | 执行 Python 代码块 | 需要 Python 逻辑处理 |
| `script:` | 执行外部脚本 | 复杂分析（R 脚本、Python 脚本） |
| `notebook:` | 执行 Jupyter notebook | 交互式分析 |

**本项目统一模板** — 每个 rule 包含：`conda` + `container` + `log` + `benchmark` + `params` + `shell`

```python
# hisat2.smk:8-28 — Rule 完整示例
rule hisat2_index:
    input:
        fasta=config["reference"]["genome_fasta"],
    output:
        protected(multiext("workflow/data/ref/hisat2_index/genome",
            ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
            ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")),
    cache: True,
    params:
        prefix=config["reference"]["hisat2_index_prefix"],
    conda: "../envs/hisat2.yaml"
    container: CONTAINER_URI
    log:     os.path.join(LOG_ROOT, "hisat2_index.log")
    benchmark: os.path.join(BENCH_ROOT, "hisat2_index.txt")
    threads: 4
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {params.prefix}) 2>/dev/null || true
        hisat2-build -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """
```
> 🔗 [hisat2.smk:8-28](workflow/rules/hisat2.smk#L8-L28) — `shell:` 完整示例（conda+container+log+benchmark+params+shell）

### R 脚本调用模式

```python
# deseq2.smk:9-45 — shell 调用外部 R 脚本的完整模式
rule deseq2:
    input:
        counts=os.path.join(OUT_FEATURECOUNTS, "featurecounts.txt"),
        gtf=config["reference"]["annotation_gtf"],
    output:
        *get_deseq2_outputs(config),
    conda:
        "../envs/deseq2.yaml"
    container:
        CONTAINER_URI
    log:
        os.path.join(LOG_ROOT, "deseq2.log"),
    params:
        groups_json=get_deseq2_groups_json(config),     # JSON 序列化传给 R
        contrasts_json=get_deseq2_contrasts_json(config),
        outdir=config["paths"]["outputs"]["deseq2"],
        padj=config["deseq2"]["padj_threshold"],
        log2fc=config["deseq2"]["log2fc_threshold"],
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        Rscript workflow/scripts/deseq2.R \
            --counts {input.counts} \
            --groups '{params.groups_json}' \
            --contrasts '{params.contrasts_json}' \
            --outdir {params.outdir} \
            --gtf {input.gtf} \
            --padj {params.padj} \
            --log2fc {params.log2fc} \
            > {log} 2>&1
        """
```
> 🔗 [deseq2.smk:9-45](workflow/rules/deseq2.smk#L9-L45) — shell 调用外部 R 脚本的完整模式

注意 JSON 参数用**单引号**包裹 `'{params.groups_json}'`，防止 shell 解析 JSON 内部的 `"` 双引号。

### PE 模式下的 Shell 管道

```bash
# hisat2.smk:78-88 — PE 模式 shell 管道: hisat2 | samtools sort
set -euo pipefail
mkdir -p {config[paths][outputs][hisat2]}
hisat2 -x {params.index_prefix} \
    -1 {input.r1} -2 {input.r2} \
    -p {threads} {params.extra} \
    2> {log} \
    | samtools sort -@ {threads} -o {output.bam} -
samtools index {output.bam}
```
> 🔗 [hisat2.smk:62-88](workflow/rules/hisat2.smk#L62-L88) — PE 模式 shell 管道

### ⭐ 最佳实践

- **`set -euo pipefail` 必须放在 shell 块首行**：确保管道中任一命令失败都能被捕获
- **output 声明要完整**：漏掉一个 output 文件，Snakemake 不知道它被产出，下次还会重跑
- **JSON 参数用单引号包裹**：防止 shell 解析 JSON 双引号造成参数截断
- **GTF 放在 `input:` 而非 `params:`**：确保 GTF 更新后 Snakemake 自动重跑下游

### ⚠ 常见错误

- **output 使用 `expand()`**：`expand()` 会展开为所有可能组合，output 应该是**本 rule 一次运行**产生的文件
- **忘记 `set -euo pipefail`**：管道 `hisat2 | samtools sort` 中 hisat2 失败时，没有 `pipefail` 的话 shell 只检查最后一个命令（samtools）的退出码
- **在 shell 中引用不存在的 input key**：PE/SE 切换时 input key 名称不同，要在 rule 级别做好切换

### 💡 面试参考 → 附录 D-3, D-4

---

## 03. 开发工作流与调试

### 开发工作流

```
Dry-run 校验（-n）→ 小范围测试（单样本）→ 全量运行 → 查看报告 → 优化
```

### 调试命令

#### `-n` / `-np`：Dry-run 预览

```bash
# 预览 DAG 执行计划（不实际运行）
$ snakemake -s workflow/Snakefile -np

# 仅预览，不打印 shell 命令
$ snakemake -s workflow/Snakefile -n
```

输出示例：
```
Job stats:
job              count
-----            -----
fastqc_raw       6
trimmomatic      6
hisat2_align     6
featurecounts    1
...
total            26
```

#### `-p`：打印实际 shell 命令

```bash
$ snakemake -s workflow/Snakefile -p
```

每个 job 执行前打印完整的 shell 命令，方便调试参数传递。

#### `--summary`：文件状态表

```bash
# 查看所有 output 文件的状态
$ snakemake -s workflow/Snakefile --summary

# 结合 --forcerun 查看重跑后哪些文件会更新
$ snakemake -s workflow/Snakefile --summary --forcerun trimmomatic
```

输出列：`output_file`, `date`, `rule`, `version`, `status`, `plan`

Status 含义：
- `updated` — 文件存在且是最新的
- `to-run` — 需要运行（input 更新或 output 缺失）
- `missing` — output 缺失，将触发运行
- `incomplete` — 文件存在但不完整（中断产物）

#### `--dag`：DAG 可视化

```bash
$ snakemake -s workflow/Snakefile --dag | dot -Tpdf > dag.pdf
```

#### `--forcerun`：强制重跑

```bash
# 强制重跑 trimmomatic 及其所有下游
$ snakemake -s workflow/Snakefile --forcerun trimmomatic

# 强制重跑特定样本的 trimmomatic（通过目标文件）
$ snakemake -s workflow/Snakefile --forcerun trimmomatic \
    runs/20260701_WT_KO/results/02_trimmomatic/WT_1_R1.trimmed.fastq.gz
```

#### `--rerun-incomplete`：中断恢复

```bash
$ snakemake -s workflow/Snakefile --rerun-incomplete -j 8
```

Snakemake 自动检测未完整写出的 output 文件并重新运行对应 rule。

#### `--report`：生成 HTML 分析报告

```bash
$ snakemake -s workflow/Snakefile --report report.html
```

#### 调度策略选择

```bash
# ILP 调度器（默认，适合复杂 DAG）
$ snakemake -s workflow/Snakefile --scheduler ilp

# Greedy 调度器（适合简单线性 DAG）
$ snakemake -s workflow/Snakefile --scheduler greedy
```

#### `--touch`：标记文件为已完成

```bash
# 标记 hisat2_index 的所有 output 为已完成（不实际运行）
$ snakemake -s workflow/Snakefile --touch hisat2_index
```

#### `--lint`：语法检查

```bash
$ snakemake -s workflow/Snakefile --lint
```

### ⭐ 最佳实践

- **每次改完 Snakefile 先 `-np`**：确认 DAG 结构符合预期
- **`--summary` 是诊断"为什么不跑"的第一工具**：快速定位 output 缺失或 input 更新
- **中断后用 `--rerun-incomplete` 而非重头跑**：避免重复已完成的计算
- **`-p` + `-n` 组合排查参数**：确认 shell 中变量替换结果正确

### ⚠ 常见错误

- **`-n` 不打印 shell 命令**：需要 `-np` 组合才能看到完整的 shell 内容
- **`--forcerun` 只强制指定 rule 重跑，不强制其上游**：上游若已满足则跳过
- **`--touch` 不更新文件内容**：只是修改时间戳，适用于重建索引等免跑场景

---

# 第二部分：核心概念

> 目标：掌握 Snakemake 的核心机制 — wildcard、expand、config、环境、日志、Python 集成

## 04. Wildcards

### Wildcard 原理

`{sample}` 是 Snakemake 的**通配符**（wildcard），在运行时自动匹配文件路径中的变量部分。

```python
# hisat2.smk:39-41 — {sample} 在 output 路径中
output:
    bam=os.path.join(OUT_HISAT2, "{sample}.sorted.bam"),
    bai=os.path.join(OUT_HISAT2, "{sample}.sorted.bam.bai"),
```
> 🔗 [hisat2.smk:39-41](workflow/rules/hisat2.smk#L39-L41)

当 Snakemake 需要 `runs/.../03_hisat2/WT_1.sorted.bam` 时，自动将 `{sample}` 匹配为 `WT_1`。这个值同时传递给同 rule 的 `input:`、`params:`、`log:`、`benchmark:` 中的 `{sample}`。

### 多 Wildcard

```python
# 多 wildcard 示例
rule process_lane:
    output:
        "results/{sample}_{lane}.bam",
```

### `wildcard_constraints`：限制通配符取值范围

```python
rule extract_chr:
    output:
        "results/{chromosome}.bam",
    wildcard_constraints:
        chromosome = "chr[0-9]+|chrX|chrY|chrM",
```

### `glob_wildcards()`：运行时动态发现 wildcard 值

```python
# 无需提前知道样本列表，从已有文件自动推断
samples, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")
```

### input 函数访问 wildcard

```python
# fastqc.smk:13 — lambda 接收 wildcards 参数
input:
    lambda wildcards: get_read_input_list(wildcards, config),
```
> 🔗 [fastqc.smk:13](workflow/rules/fastqc.smk#L13)

input 函数签名为 `(wildcards) -> str|list|dict`，每次为新 `{sample}` 调用时，`wildcards.sample` 即为当前样本名。

### `unpack()` — 将 dict key 展开为独立 input

```python
# hisat2.smk:37-38 — unpack() 将 dict 展开为 input.read 等
input:
    unpack(lambda wildcards: get_hisat2_input(wildcards, config)),
    # get_hisat2_input 返回 {"r1": ..., "r2": ...} (PE) 或 {"read": ...} (SE)
    # unpack 后 shell 中通过 {input.r1}, {input.r2} 或 {input.read} 引用
```
> 🔗 [hisat2.smk:37-38](workflow/rules/hisat2.smk#L37-L38)

### Wildcard 的取值范围

```python
# Snakefile:79 — SAMPLES 列表定义了 {sample} 的所有可能值
SAMPLES = get_samples(config)  # ["WT_1", "WT_2", "KO_1", ...]
```
> 🔗 [Snakefile:79](workflow/Snakefile#L79)

### ⭐ 最佳实践

- **output 和 input 中 wildcard 名称必须一致**：output 中 `{sample}` 对应 input 中 `{sample}`
- **用 `glob_wildcards()` 发现样本**：避免手动维护样本列表与文件系统不同步
- **`unpack()` + λ 函数**是 PE/SE 切换的关键模式：不同模式下 input key 不同，`unpack()` 自动处理

### ⚠ 常见错误

- **output 用 `{sample}` 但 input 用了 `{s}`**：名称不匹配会导致 Snakemake 无法推导 wildcard 值
- **在 input 函数外使用 `wildcards` 变量**：wildcards 只在 input/output 的函数闭包内可用
- **未定义 wildcard 取值范围**：Snakemake 从 `rule all` 的目标文件名中自动推断，但如果规则产出不在 `rule all` 中，需要确保目标由其他规则"请求"

### 💡 面试参考 → 附录 D-5, D-6

---

## 05. expand()

### expand 原理：笛卡尔积生成文件列表

`expand()` 是 Snakemake 内置函数，生成文件路径的**笛卡尔积**：

```python
expand("results/{sample}_{read}.html", sample=["WT_1", "KO_1"], read=["R1", "R2"])
# → ["results/WT_1_R1.html", "results/WT_1_R2.html",
#    "results/KO_1_R1.html", "results/KO_1_R2.html"]
```

### expand 用于 input 收集已有文件（不是 output 声明）

**关键原则**：`expand()` 用于 `input:` 收集已有文件，**不能**用于 `output:` 声明产出。

```python
# multiqc.smk:10-20 — expand() + 条件表达式适配 PE/SE 双模式
input:
    fastqc_raw=expand(
        [os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.html"),
         os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.html"),
         os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.zip"),
         os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.zip")],
        sample=SAMPLES,
    ),
```
> 🔗 [multiqc.smk:10-20](workflow/rules/multiqc.smk#L10-L20) — `expand()` 收集所有样本的 FastQC 输出

```python
# multiqc.smk:32 — expand() 收集所有样本的 HISAT2 日志
hisat2_logs=expand(f"{LOG_ROOT}/hisat2_align/{{sample}}.log", sample=SAMPLES),
```
> 🔗 [multiqc.smk:32](workflow/rules/multiqc.smk#L32)

注意：expand 模板中的 `{sample}` 需要双花括号 `{{sample}}` 转义，区分于 expand 的参数 `sample=SAMPLES`。

### 列表推导式等效写法

```python
# featurecounts.smk:11 — 列表推导式（与 expand 等效）
bams=[f"{OUT_HISAT2}/{s}.sorted.bam" for s in SAMPLES],
```
> 🔗 [featurecounts.smk:11](workflow/rules/featurecounts.smk#L11)

`expand()` 和列表推导式的选择：
- `expand()` — 语义更清晰，适合多维展开（笛卡尔积）
- 列表推导式 — 更灵活，适合单维或需要条件判断的场景

### 多维 expand 与 zip 模式

```python
# 笛卡尔积：sample × read → 所有组合
expand("data/{sample}_{read}.fq", sample=SAMPLES, read=["R1", "R2"])

# zip 模式：一一对应，不产生笛卡尔积
expand("data/{sample}_{read}.fq", zip, sample=SAMPLES, read=["R1", "R2"])
# SAMPLES=["A","B"], read=["R1","R2"] → ["data/A_R1.fq", "data/B_R2.fq"]
```

### ⭐ 最佳实践

- **expand 用于 input，不用于 output**：output 应声明本 rule 一次运行产出的文件
- **expand 中的 wildcard 用双花括号转义**：`{sample}` → `{{sample}}`
- **expand 理解成"我要收集这些文件作为 input"**，不是"我要生成这些文件"

### ⚠ 常见错误

- **在 output 中使用 expand**：output 需要用 `{wildcard}` 声明单个产出，而非笛卡尔积展开
- **expand 笛卡尔积陷阱**：`expand("data/{sample}_{read}.html", sample=SAMPLES, read=READS)` 会产生 N×M 个路径，确认这是你想要的行为
- **单花括号导致与 expand 参数混淆**：模板中的 `{sample}` 被当作 expand 的参数名而非 wildcard

### 💡 面试参考 → 附录 D-7, D-8

---

## 06. Config 管理（含 params / resources）

### `config.yaml` 结构与 `configfile:` 指令

```python
# Snakefile:70 — 加载配置文件
configfile: "config/config.yaml"
```
> 🔗 [Snakefile:70](workflow/Snakefile#L70)

`configfile:` 指令将 YAML 文件加载为 Python dict，后续代码通过 `config["key"]` 访问。

### 项目路径配置：批次隔离 `{batch}` 占位符

```python
# Snakefile:73-76 — 批次路径格式化
_batch = config["paths"]["batch"]
config["paths"]["outputs"] = {k: v.format(batch=_batch) for k, v in config["paths"]["outputs"].items()}
config["paths"]["log_root"] = config["paths"]["logs"].format(batch=_batch)
config["paths"]["bench_root"] = config["paths"]["benchmarks"].format(batch=_batch)
```
> 🔗 [Snakefile:73-76](workflow/Snakefile#L73-L76)

```yaml
# config.yaml:11-23 — batch 占位符在 output 模板中
paths:
  batch: "20260701_WT_KO"
  outputs:
    fastqc_raw: "runs/{batch}/results/01_fastqc_raw"
    trimmomatic: "runs/{batch}/results/02_trimmomatic"
    # ...
```
> 🔗 [config.yaml:11-23](config/config.yaml#L11-L23)

### 使用场景表

| 使用场景 | 在哪 | 示例 |
|---------|------|------|
| 定义默认值 | config.yaml | `batch: "20260701_WT_KO"` |
| 控制流程分支 | .smk 规则 | `if config["read_pattern"]["mode"] == "single"` |
| 透传给 shell | params | `params: extra=config["hisat2"]["extra"]` |
| JSON 序列化给 R | params | `params: groups_json=get_deseq2_groups_json(config)` |
| shell 中直接引用 | shell 块 | `{config[paths][outputs][hisat2]}` |
| 资源分级分配 | profile | `set-resources: { hisat2_index: {mem_mb: 32000} }` |
| 启动前校验 | common.py | `validate_config(config)` |

### `params` — 从 config 透传参数

```python
# deseq2.smk:23-29 — params 从 config 取值，序列化为 JSON 传给 R 脚本
params:
    groups_json=get_deseq2_groups_json(config),
    contrasts_json=get_deseq2_contrasts_json(config),
    outdir=config["paths"]["outputs"]["deseq2"],
    padj=config["deseq2"]["padj_threshold"],
    log2fc=config["deseq2"]["log2fc_threshold"],
    top_n=config["deseq2"]["top_n_genes"],
```
> 🔗 [deseq2.smk:23-29](workflow/rules/deseq2.smk#L23-L29)

### `resources` — 资源限制 + `set-resources` 分层分配

```yaml
# slurm-apptainer/config.yaml:30-59 — default-resources + set-resources 分层资源分配
default-resources:
  slurm_partition: "compute"
  runtime: 240
  mem_mb: 16000
  cpus_per_task: 4

set-resources:
  hisat2_index:
    runtime: 120
    mem_mb: 32000
    cpus_per_task: 8
  hisat2_align:
    runtime: 360
    mem_mb: 32000
    cpus_per_task: 8
```
> 🔗 [slurm-apptainer/config.yaml:30-59](profile/slurm-apptainer/config.yaml#L30-L59)

### `threads` — Snakemake 内置，自动传给调度器

```python
# hisat2.smk:22 — threads 在规则中直接引用 {threads}
threads: 4
# shell 中: hisat2-build -p {threads} ...
```

### ⭐ 最佳实践

- **唯一配置入口**：所有参数集中在 `config.yaml`，不要在 .smk 文件中硬编码
- **params 中转**：config 值通过 params 传给 shell，不要在 shell 中直接写 `config["xxx"]`（除非是路径类）
- **启动前校验**：`validate_config()` 检查关键字段，避免跑到一半才发现拼写错误
- **批次隔离**：改 `paths.batch` 一个值，所有产出自动切换到新目录

### ⚠ 常见错误

- **`config[key]` 在 Python 层直接访问，shell 层用 `{config[key]}`**：两者语法不同
- **params 和 resources 混用**：params 是参数值，resources 是资源限制（内存/时间/CPU）
- **修改 config.yaml 后不检查**：用 `--lint` + `-np` 验证语法

### 💡 面试参考 → 附录 D-9, D-10

---

## 07. Conda 环境管理

### 通道顺序规范

Bioconda 官方推荐通道顺序：**conda-forge → bioconda → defaults**。避免 conda-forge 的新版本被 bioconda 旧版本覆盖。

```yaml
# envs/hisat2.yaml:1-9 — Conda 环境定义
name: hisat2_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - hisat2=2.2.1
  - samtools=1.18
```
> 🔗 [envs/hisat2.yaml:1-9](workflow/envs/hisat2.yaml)

### 版本固定（Pin 具体版本号）

精确指定版本号（`hisat2=2.2.1` 而非 `hisat2`）确保可复现。本项目所有 conda 环境均固定版本：

| 工具 | 版本 | Conda 环境 |
|------|------|-----------|
| FastQC | 0.12.1 | envs/fastqc.yaml |
| Trimmomatic | 0.39 | envs/trimmomatic.yaml |
| HISAT2 | 2.2.1 | envs/hisat2.yaml |
| Samtools | 1.18 | envs/hisat2.yaml |
| featureCounts | 2.0.6 | envs/featurecounts.yaml |
| MultiQC | 1.21 | envs/multiqc.yaml |
| R | 4.3.2 | envs/deseq2.yaml |
| DESeq2 | 1.42.0 | envs/deseq2.yaml |
| clusterProfiler | 4.10.0 | envs/clusterprofiler.yaml |

### 每个模块独立环境 vs 统合环境

本项目采用**每模块独立环境**策略：每个 rule 有专属 `.yaml`。这样做的好处是环境更小（安装快、冲突少），但需要管理多个文件。

容器化生产环境则使用**统合镜像**（`containers/Dockerfile`），包含所有工具，适合一次性部署。

### `conda:` 指令与路径规则

```python
# hisat2.smk:18 — conda 路径相对于 .smk 文件所在目录
conda: "../envs/hisat2.yaml"
```
> 🔗 [hisat2.smk:18](workflow/rules/hisat2.smk#L18)

`.smk` 文件在 `workflow/rules/`，conda 环境在 `workflow/envs/`，所以相对路径为 `../envs/xxx.yaml`。

### conda 与 container 互斥

```bash
# 开发模式 — 使用 Conda
$ snakemake -s workflow/Snakefile --software-deployment-method conda -j 8

# 生产模式 — 使用容器（镜像已包含全部工具，无需 conda）
$ snakemake -s workflow/Snakefile --profile profile/docker
```

`--software-deployment-method conda` 时不激活 container；容器模式设置 `software-deployment-method: [apptainer]`（镜像已包含全部工具）。

### ⭐ 最佳实践

- **通道顺序固定**：conda-forge → bioconda → defaults
- **版本精确固定**：Pin 到具体版本号（`=1.2.3`），不依赖 latest
- **独立环境 vs 统合镜像分工**：开发用独立 conda 环境（快速迭代），生产用统合容器镜像（一次部署）
- **环境文件与 rule 放在同一仓库**：相对路径引用，clone 即可用

### ⚠ 常见错误

- **通道顺序错误**：defaults → bioconda → conda-forge 可能导致依赖解析失败
- **相对路径写错**：conda 路径相对于 `.smk` 文件，不是相对于 Snakefile
- **开发环境与生产环境不一致**：conda 环境和容器镜像的工具版本需要保持同步

### 💡 面试参考 → 附录 D-11

---

## 08. 日志与性能监控（Log + Benchmark）

### Log 管理

Snakemake 的 `log:` directive 自动为每个 job 创建独立的日志文件，并自动创建父目录。

```python
# hisat2.smk:44 — {sample} wildcard 在 log 路径中
log: os.path.join(LOG_ROOT, "hisat2_align/{sample}.log"),
```
> 🔗 [hisat2.smk:44](workflow/rules/hisat2.smk#L44)

**重要**：`mkdir -p` 只用于 output 目录，**不用于 log 目录**——Snakemake 自动为 `log:` 创建父目录。

```python
# Snakefile:97-98 — LOG_ROOT 和 BENCH_ROOT 常量定义
LOG_ROOT   = config["paths"]["log_root"]
BENCH_ROOT = config["paths"]["bench_root"]
```
> 🔗 [Snakefile:97-98](workflow/Snakefile#L97-L98)

### R 脚本日志规范

三份 R 脚本统一使用 `[模块名]` 前缀消息格式：

```r
# deseq2.R:59-66 — R 脚本 [DESeq2] 模块前缀日志
message("[DESeq2] ========================================")
message("[DESeq2] DESeq2 差异表达分析")
message("[DESeq2] ========================================")
message("[DESeq2] Counts file: ", counts_file)
message("[DESeq2] Output dir:  ", outdir)
message("[DESeq2] padj cutoff: ", padj_cut)
message("[DESeq2] log2FC cutoff: ", log2fc_cut)
message("[DESeq2]")
```
> 🔗 [deseq2.R:59-66](workflow/scripts/deseq2.R#L59-L66)

统一前缀便于 `grep` 从 Snakemake 日志中过滤特定模块输出：
```bash
$ grep "\[DESeq2\]" runs/20260701_WT_KO/logs/deseq2.log
```

### Benchmark

`benchmark:` directive 记录每次运行的性能数据。

```python
# featurecounts.smk:22-23 — benchmark 定义
benchmark:
    os.path.join(BENCH_ROOT, "featurecounts.txt"),
```
> 🔗 [featurecounts.smk:22-23](workflow/rules/featurecounts.smk#L22-L23)

Benchmark 输出字段：

| 字段 | 含义 | 用途 |
|------|------|------|
| `s` | CPU 时间（秒） | 计算密集型程度 |
| `wall_clock` | 实际耗时（秒） | 用户感知的等待时间 |
| `max_rss` | 最大物理内存（MB） | 内存资源需求评估 |
| `max_vms` | 最大虚拟内存（MB） | 内存地址空间使用 |
| `io_in` / `io_out` | 磁盘读写量（MB） | IO 密集型程度 |
| `mean_load` | 平均 CPU 负载 | 并行效率评估 |
| `cpu_time` | 用户态+内核态 CPU 时间 | 与 threads 对比评估并行效率 |

### 如何用 benchmark 定位慢任务

```bash
# 找出耗时最长的 rule
$ grep -l "" runs/*/benchmarks/*.txt | xargs grep "wall_clock" | sort -t: -k2 -rn | head -10

# 找内存消耗最大的 job
$ grep "max_rss" runs/*/benchmarks/*.txt | sort -t: -k2 -rn | head -5
```

### ⭐ 最佳实践

- **每个 rule 都有独立的 log + benchmark**：这是本项目模板的标准
- **log 路径包含 `{sample}`**：每个样本的日志独立，便于并行排查
- **R 脚本日志用模块前缀**：`[DESeq2]` / `[clusterProfiler]` 方便过滤
- **定期查看 benchmark**：识别资源瓶颈，优化 `set-resources` 配置

### ⚠ 常见错误

- **在 shell 中 `mkdir -p` log 目录**：Snakemake 自动创建，多此一举还容易路径不一致
- **benchmark 路径写死**：应使用 `BENCH_ROOT` 常量，保持与项目路径体系一致
- **R 脚本用 `print()` 而非 `message()`**：`message()` 输出到 stderr 被 `> {log} 2>&1` 捕获

---

## 09. Python 函数与 Snakemake 集成

### 函数在 Snakefile 中的角色

Python 函数在三个位置被 Snakemake 调用：`input:` / `output:` / `params:`。本项目将所有辅助函数提取到 `common.py`，实现**单一职责 + 模块化设计**。

### `common.py` 函数职责速查

| 使用场景 | 函数 | 被调用位置 |
|---------|------|-----------|
| 获取样本列表 | `get_samples(config)` | [Snakefile:79](workflow/Snakefile#L79) |
| PE/SE 感知 input | `get_read_inputs(wildcards, config)` | [trimmomatic.smk:12](workflow/rules/trimmomatic.smk#L12) |
| 动态 output 列表 | `get_deseq2_outputs(config)` | [deseq2.smk:14](workflow/rules/deseq2.smk#L14) |
| rule all 目标汇总 | `get_all_pipeline_targets(config, samples)` | [Snakefile:123](workflow/Snakefile#L123) |
| 启动前校验 | `validate_config(config)` | [Snakefile:82](workflow/Snakefile#L82) |
| JSON 序列化 | `to_json_str(obj)` | [common.py:80](workflow/scripts/common.py#L80) |

### input 函数 λ 闭包传送 config

```python
# common.py:88-95 — get_read_input_list() 模式感知的 λ 函数
def get_read_input_list(wildcards, config: dict) -> List[str]:
    mode = _get_read_mode(config)
    raw_dir = config["paths"]["raw_dir"]
    pattern = config["read_pattern"]
    if mode == "single":
        return [os.path.join(raw_dir, f"{wildcards.sample}{pattern['single_suffix']}")]
    r1, r2 = get_read_pairs(wildcards.sample, raw_dir, pattern["r1_suffix"], pattern["r2_suffix"])
    return [r1, r2]
```
> 🔗 [common.py:88-95](workflow/scripts/common.py#L88-L95)

在 .smk 中通过 λ 闭包传递 config：

```python
# fastqc.smk:13 — lambda 闭包传递 config
input:
    lambda wildcards: get_read_input_list(wildcards, config),
```

### `*get_xxx_outputs()` Star unpacking 动态 output

```python
# common.py:240-254 — get_deseq2_outputs() 按对比动态展开输出列表
def get_deseq2_outputs(config: dict) -> List[str]:
    deseq2_cfg = config.get("deseq2", {})
    contrasts = deseq2_cfg.get("contrasts", [])
    outdir = _out(config, "deseq2")
    outputs = [
        f"{outdir}/PCA_plot.pdf", f"{outdir}/sample_distance_heatmap.pdf",
        f"{outdir}/DEG_heatmap.pdf",
    ]
    for ct in contrasts:
        n = ct["name"]
        outputs.extend([
            f"{outdir}/{n}_all_results.csv",  f"{outdir}/{n}_significant.csv",
            f"{outdir}/{n}_MA_plot.pdf",      f"{outdir}/{n}_volcano_plot.pdf",
        ])
    return outputs
```
> 🔗 [common.py:240-254](workflow/scripts/common.py#L240-L254)

```python
# deseq2.smk:14 — *get_deseq2_outputs(config) Star unpacking
output:
    *get_deseq2_outputs(config),
```
> 🔗 [deseq2.smk:14](workflow/rules/deseq2.smk#L14)

`*` 将 Python 列表展开为多个 output 项，Snakemake 据此跟踪每个文件。

### `unpack()` 动态 input key

```python
# hisat2.smk:37-38 — unpack() 动态 input key
input:
    unpack(lambda wildcards: get_hisat2_input(wildcards, config)),
```
> 🔗 [hisat2.smk:37-38](workflow/rules/hisat2.smk#L37-L38)

`get_hisat2_input()` 返回 `{"r1": ..., "r2": ...}` (PE) 或 `{"read": ...}` (SE)，`unpack()` 将 dict key 展开为 shell 中的 `{input.r1}`、`{input.r2}` 或 `{input.read}`。

### 启动前校验

```python
# common.py:281-345 — validate_config() 6 类配置校验
def validate_config(config: dict) -> None:
    errors: List[str] = []

    # 1. read_pattern.mode 必须是 paired 或 single
    mode = config.get("read_pattern", {}).get("mode", "paired")
    if mode not in ("paired", "single"):
        errors.append(f"read_pattern.mode 必须为 'paired' 或 'single'，当前值: '{mode}'")

    # 2. 模式对应后缀必须存在
    # 3. 参考文件存在性
    # 4. deseq2 contrasts 结构
    # 5. clusterprofiler 段
    # 6. 样本表文件存在性

    if errors:
        print("ERROR: config.yaml 配置校验失败：", file=sys.stderr)
        for e in errors:
            print(f"  ✗ {e}", file=sys.stderr)
        raise SystemExit(1)
```
> 🔗 [common.py:281-345](workflow/scripts/common.py#L281-L345) — `validate_config()` 6 类配置校验

### ⭐ 最佳实践

- **函数单一职责**：一个函数只做一件事，`common.py` 中 16 个函数各司其职
- **input 函数用 λ 闭包传递 config**：因为函数在 `common.py` 模块中无法访问 Snakemake 全局变量
- **动态 output 用 `*list` unpacking**：新增对比后 output 自动扩展，无需手动修改 rule
- **启动前校验是必须的**：避免跑到一半才发现配置错误（如拼写错误的列名、不存在的文件路径）

### ⚠ 常见错误

- **input 函数内的逻辑太重**：input 函数每次 job 都被调用，避免大量 IO 或计算
- **在 common.py 中直接访问 config 全局变量**：config 通过参数传入，保持函数纯净化
- **Star unpacking 返回空列表**：如果没有对比定义，`*[]` 展开为空，rule 变成无 output（Snakemake 报错）

---

# 第三部分：实战进阶

> 目标：多样本调度、文件标记、冲突处理、checkpoint、HPC、生命周期

## 10. 多样本调度

### `{sample}` wildcard 自动并行化

Snakemake 为每个 `{sample}` 值生成独立的 job。同一层级无依赖的 job 自动并行：

```
hisat2_align (WT_1)  ──┐
hisat2_align (WT_2)  ──┼──→ featurecounts
hisat2_align (KO_1)  ──┤
hisat2_align (KO_2)  ──┘
```

4 个样本的 `hisat2_align` 无相互依赖，`-j 4` 即可完全并行。

### PE/SE 在 rule 定义层通过 `if/else` 切换

```python
# fastqc.smk:10-53 — PE/SE 在 rule 级切换（两个完整 rule 定义）
if config["read_pattern"]["mode"] == "single":
    rule fastqc_raw:
        input:
            lambda wildcards: get_read_input_list(wildcards, config),
        output:
            html=os.path.join(OUT_FASTQC_RAW, "{sample}_fastqc.html"),
            zip=os.path.join(OUT_FASTQC_RAW,  "{sample}_fastqc.zip"),
        # ... single-end specific config ...
else:
    rule fastqc_raw:
        input:
            lambda wildcards: get_read_input_list(wildcards, config),
        output:
            html_r1=os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.html"),
            html_r2=os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.html"),
            zip_r1=os.path.join(OUT_FASTQC_RAW,  "{sample}_R1_fastqc.zip"),
            zip_r2=os.path.join(OUT_FASTQC_RAW,  "{sample}_R2_fastqc.zip"),
        # ... paired-end specific config ...
```
> 🔗 [fastqc.smk:10-53](workflow/rules/fastqc.smk#L10-L53) — PE/SE 在 rule 级切换

```bash
# hisat2.smk:51-59 — SE: -U 单端输入
hisat2 -x {params.index_prefix} \
    -U {input.read} \
    -p {threads} {params.extra} \
    2> {log} \
    | samtools sort -@ {threads} -o {output.bam} -

# hisat2.smk:80-87 — PE: -1/-2 双端输入
hisat2 -x {params.index_prefix} \
    -1 {input.r1} -2 {input.r2} \
    -p {threads} {params.extra} \
    2> {log} \
    | samtools sort -@ {threads} -o {output.bam} -
```
> 🔗 [hisat2.smk:51-59](workflow/rules/hisat2.smk#L51-L59) — SE 模式
> 🔗 [hisat2.smk:80-87](workflow/rules/hisat2.smk#L80-L87) — PE 模式

关键设计决策：**在 rule 级切换而非 shell 内部 `if/else`**。原因：
1. Snakemake 在 DAG 构建阶段需要知道 input/output 的完整列表
2. shell 中的条件分支不会改变 Snakemake 看到的 input/output 声明
3. PE/SE 的 input key 名称不同（`r1/r2` vs `read`），shell 模板要求所有引用键必须存在

PE/SE 切换只需改一个值：
```yaml
# config.yaml:46-47 — 一个值控制全局
read_pattern:
  mode: "paired"   # ← 改为 "single" 即可切换
```
> 🔗 [config.yaml:46-47](config/config.yaml#L46-L47)

### `rule all` 动态生成目标列表

```python
# Snakefile:114-123 — rule all 调用 Python 函数动态生成目标
rule all:
    input:
        get_all_pipeline_targets(config, SAMPLES),
```
> 🔗 [Snakefile:114-123](workflow/Snakefile#L114-L123)

`get_all_pipeline_targets()` 根据 PE/SE 模式生成不同的目标文件列表，确保切换模式后 `rule all` 自动适配。

### shadow 执行模式

```python
rule some_rule:
    shadow: "full"   # 在临时目录运行，防止并行写冲突
```

适用于多个 job 写入同一目录时可能冲突的场景。

### ⭐ 最佳实践

- **PE/SE 在 rule 级别的 `if/else` 切换**：而非 shell 内部条件判断
- **`rule all` 动态生成目标**：切换模式后无需手动修改
- **`-j` 参数合理设置**：等于/略小于可用 CPU 核数，避免资源竞争

### ⚠ 常见错误

- **PE/SE 切换后不改 config**：只要改 `mode: "single"` 即可，修改 Snakefile 是错的
- **在一条 rule 的 shell 里用 if 判断 PE/SE**：input key 在 DAG 构建时已固定，shell 中无法动态切换 key
- **`-j` 设太大导致 IO 瓶颈**：并行度受限于磁盘 IO 而非 CPU

---

## 11. temp() / protected() / directory()

### `temp()` — 中间文件，下游消费后自动删除

```python
# trimmomatic.smk:14 — temp() 标记修剪后 FASTQ 为临时文件
output:
    read=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}.trimmed.fastq.gz")),
```
> 🔗 [trimmomatic.smk:14](workflow/rules/trimmomatic.smk#L14)

PE 模式：
```python
# trimmomatic.smk:45-48 — PE 模式 4 个 temp() 输出
output:
    r1=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R1.trimmed.fastq.gz")),
    r2=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R2.trimmed.fastq.gz")),
    r1_unpaired=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R1.unpaired.fastq.gz")),
    r2_unpaired=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R2.unpaired.fastq.gz")),
```
> 🔗 [trimmomatic.smk:45-48](workflow/rules/trimmomatic.smk#L45-L48)

**删除时机**：当所有以该文件为 input 的下游 rule 都成功完成后，Snakemake 自动删除 `temp()` 文件。

**为什么 Trimmomatic 输出用 temp()？** 修剪后的 FASTQ 体积大（与原始数据相当），比对完成后不再需要，自动清理可节省大量磁盘空间。

### `protected()` — 写保护

```python
# hisat2.smk:12-14 — protected() 保护索引文件不被意外覆盖
output:
    protected(multiext("workflow/data/ref/hisat2_index/genome",
        ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
        ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")),
```
> 🔗 [hisat2.smk:12-14](workflow/rules/hisat2.smk#L12-L14)

`protected()` 对 output 文件加写保护（`chmod a-w`）。要强制覆盖需 `--forceall --protected` 参数。

**为什么索引用 protected()？** 基因组索引构建耗时（可能数小时），意外覆盖代价大。

### `multiext()` — 多扩展名输出声明

```python
multiext("path/to/genome", ".1.ht2", ".2.ht2", ..., ".8.ht2")
# 展开为:
# ["path/to/genome.1.ht2", "path/to/genome.2.ht2", ..., "path/to/genome.8.ht2"]
```

HISAT2 索引由 8 个 `.ht2` 文件组成，`multiext()` 显式声明全部文件。替代旧式的 `touch()` 哨兵文件模式——任意一个 `.ht2` 文件缺失，Snakemake 都会检测到并重跑。

### `cache: True` — 内容寻址跨批次缓存

```python
# hisat2.smk:15 — 内容寻址缓存
cache: True,
```
> 🔗 [hisat2.smk:15](workflow/rules/hisat2.smk#L15)

与时间戳缓存不同，`cache: True` 基于文件**内容哈希**判断是否需要重跑。同一个基因组不同批次间可跨批次复用索引，无需重新构建。

### `directory()` — 目录级输出

```python
rule some_tool:
    output:
        directory("results/my_output_dir/"),
    shell:
        "some_tool --outdir {output}"
```

Shell 需自行保证目录内容完整（Snakemake 只检查目录是否存在）。

### `ancient()` — 标记永不更新的输入

```python
input:
    ancient("data/ref/genome.fa"),
```

时间戳改变时不会触发重跑（适用于参考基因组等稳定文件）。

### ⭐ 最佳实践

- **中间大文件用 temp()**：自动清理节省磁盘
- **构建耗时的产出用 protected()**：防止误删或意外覆盖
- **多文件产出用 multiext()**：显式声明全部文件，替代 `touch()` 哨兵
- **跨批次可复用的产出用 cache: True**：内容寻址 > 时间戳

### ⚠ 常见错误

- **temp() 文件在内存中被下游 rule 依赖时不会删除**：只有所有下游全完成才删
- **protected() 文件忘记 `--forceall --protected` 参数无法覆盖**：需要同时指定两个 flag
- **`multiext()` 声明不全**：遗漏扩展名会导致部分文件缺失时 Snakemake 不感知

### 💡 面试参考 → 附录 D-12, D-13

---

## 12. ruleorder

### 概念

当多个 rule 能产生同一个 output 时，Snakemake 会报歧义错误（AmbiguousRuleException）。`ruleorder` 消除歧义。

### 语法

```python
ruleorder: specific_rule > generic_rule
```

### 典型场景

不同样本使用不同的比对流程：

```python
rule hisat2_align_special:
    output:
        "results/align/{sample}.bam",
    # 特殊样本用 STAR 比对...

rule hisat2_align:
    output:
        "results/align/{sample}.bam",
    # 通用样本用 HISAT2 比对...

ruleorder: hisat2_align_special > hisat2_align
```

> 🔗 本项目当前所有 output 仅一个 rule 产生，不需要 ruleorder。以上为假设场景示例。

### ⭐ 最佳实践

- **优先设计互斥的 output 路径**（如 `results/hisat2/{sample}.bam` vs `results/star/{sample}.bam`），自然避免歧义
- **ruleorder 是最后手段**：路径设计解决 > 参数化解决 > ruleorder 解决

### ⚠ 常见错误

- **两个 rule 产出相同路径但 rule 名称不同**：Snakemake 不关心 rule 名称，只关心 output 路径
- **ruleorder 写反**：`specific > generic` 不是 `generic > specific`

### 💡 面试参考 → 附录 D-14

---

## 13. checkpoint

### 概念

`checkpoint` 是一种特殊的 rule，其输出文件数量和名称在**运行时**才知道。普通 rule 的 output 在 DAG 构建时完全确定，checkpoint 则推迟到执行阶段。

### 与普通 rule 的区别

| 特性 | 普通 rule | checkpoint |
|------|----------|------------|
| output 确定性 | DAG 构建时已知 | 运行时才知道 |
| wildcard 来源 | 提前定义（SAMPLES） | 运行时动态产生 |
| DAG 影响 | 完整的静态 DAG | DAG 在执行中动态扩展 |

### 典型场景

- 样本拆分（demultiplexing）：一个 FASTQ 包含多个样本，运行时才知道拆出多少个
- 动态聚类：聚类数由数据决定，不同运行产出不同数量的聚类结果

### 使用方式

```python
checkpoint demultiplex:
    input:
        "data/raw/pooled.fastq.gz",
    output:
        directory("results/demux/"),
    shell:
        """
        demux_tool {input} --outdir {output}
        """

# 在 Python 代码中获取运行时结果
def get_demuxed_samples(wildcards):
    demux_output = checkpoints.demultiplex.get()  # 获取运行后的 output
    sample_ids = glob_wildcards(f"{demux_output}/{{sample}}.fastq.gz").sample
    return expand("results/align/{sample}.bam", sample=sample_ids)
```

> 🔗 本项目未使用 checkpoint。以上为假设场景示例。

### ⭐ 最佳实践

- **只在确实需要运行时才知道 output 时才用 checkpoint**：大多数场景可以用 input 函数 + `glob_wildcards()` 解决
- **checkpoint 下游规则必须通过函数访问**：不能用静态 `expand()`

### ⚠ 常见错误

- **用 checkpoint 解决可以用 input 函数解决的问题**：增加了复杂度却没有实际收益
- **忘记 `checkpoints.rule_name.get()`**：checkpoint 的 output 不能直接用 `rules.rule_name.output`

### 💡 面试参考 → 附录 D-15

---

## 14. HPC 与 Slurm

### Snakemake 9+ Executor 插件

Snakemake 9 使用插件式 executor：

```yaml
# slurm-apptainer/config.yaml:15 — executor 声明
executor: slurm
```
> 🔗 [slurm-apptainer/config.yaml:15](profile/slurm-apptainer/config.yaml#L15)

需要安装插件：
```bash
$ pip install snakemake-executor-plugin-slurm
```

### Profile 三层体系

本项目提供三套 Profile，覆盖从开发到生产的完整链路：

| Profile | 用途 | 环境 |
|---------|------|------|
| [profile/docker/](profile/docker/) | 本地开发测试 | Docker 容器 |
| [profile/apptainer/](profile/apptainer/) | 单机生产 | Apptainer 容器 |
| [profile/slurm-apptainer/](profile/slurm-apptainer/) | 集群生产 | Slurm + Apptainer |

```yaml
# profile/docker/config.yaml — Docker 本地开发配置
software-deployment-method:
  - apptainer
cores: 16
jobs: 4
latency-wait: 30
```
> 🔗 [profile/docker/config.yaml:1-21](profile/docker/config.yaml#L1-L21)

```yaml
# profile/apptainer/config.yaml — 单机 Apptainer 配置
software-deployment-method:
  - apptainer
apptainer-args: --cleanenv --bind .
apptainer-prefix: "containers/apptainer_images"
cores: 32
jobs: 8
restart-times: 1
latency-wait: 60
```
> 🔗 [profile/apptainer/config.yaml:1-39](profile/apptainer/config.yaml#L1-L39)

### `localrules:` — 标记本地规则

```python
localrules: all, create_dirs
```

标记为 `localrules` 的规则不在集群提交，直接在登录节点运行。适用于：
- `rule all`（只是声明目标）
- 目录创建等轻量任务
- 避免了为几秒的任务浪费 Slurm 调度排队时间

### `default-resources` + `set-resources` 分层资源分配

```yaml
# slurm-apptainer/config.yaml:30-59 — 分层资源分配
default-resources:          # 所有 rule 的默认值
  slurm_partition: "compute"
  runtime: 240              # 4 小时
  mem_mb: 16000
  cpus_per_task: 4

set-resources:              # 按 rule 覆盖
  hisat2_index:
    runtime: 120            # 2 小时
    mem_mb: 32000
    cpus_per_task: 8
  hisat2_align:
    runtime: 360            # 6 小时
    mem_mb: 32000
    cpus_per_task: 8
```
> 🔗 [slurm-apptainer/config.yaml:30-59](profile/slurm-apptainer/config.yaml#L30-L59)

### `jobs` 与 `latency-wait`

- `jobs: 24` — 限制同时提交到 Slurm 的最大作业数
- `latency-wait: 120` — 等待 NFS 文件系统传播文件更新的秒数（集群环境必须设置）

### Apptainer 镜像管理

```bash
# containers/build.sh — 容器构建统一入口
$ bash containers/build.sh docker         # 构建 Docker 镜像
$ bash containers/build.sh all            # Docker + 转 SIF (推荐)
$ bash containers/build.sh test           # 验证镜像内所有工具
```
> 🔗 [containers/build.sh](containers/build.sh)

### ⭐ 最佳实践

- **三层 Profile 覆盖全部场景**：本地开发（Docker）→ 单机（Apptainer）→ 集群（Slurm+Apptainer）
- **default-resources 给出合理底线**：防止某个 rule 忘记配置资源导致 OOM
- **set-resources 按 rule 的实测 benchmark 值设置**：不要猜，用 benchmark 数据
- **latency-wait 在 NFS 环境必须设置**：典型值 60-120 秒

### ⚠ 常见错误

- **`jobs` vs `cores` 混淆**：`jobs` 是 Slurm 作业数上限，`cores` 是单机并行核数
- **忘记设置 `latency-wait`**：NFS 延迟导致 output 文件已写出但下游 job 看不到
- **localrules 漏配**：`rule all` 提交到 Slurm 会浪费一个作业槽等待

### 💡 面试参考 → 附录 D-16, D-17

---

## 15. Workflow 生命周期

### 完整生命周期

```
初始化: configfile → config → SAMPLES → validate
   ↓
调度:   解析 DAG → 生成 job 队列 → 按 -j 并行度分配
   ↓
执行:   cd workdir && {shell block}（每个 job 独立运行）
   ↓
缓存:   比较 input/output 时间戳 → 决定是否运行
   ↓
完成:   所有 rule all input 满足 → 退出 0
```

```python
# Snakefile:70-83 — 初始化阶段完整流程
configfile: "config/config.yaml"                    # 1. 加载配置

_batch = config["paths"]["batch"]                   # 2. 批次占位符替换
config["paths"]["outputs"] = {k: v.format(batch=_batch)
    for k, v in config["paths"]["outputs"].items()}

SAMPLES = get_samples(config)                       # 3. 解析样本列表

validate_config(config)                             # 4. 校验配置
validate_input_files(config, SAMPLES)               # 5. 校验输入文件
```
> 🔗 [Snakefile:70-83](workflow/Snakefile#L70-L83) — 初始化阶段完整流程

### 失败处理与恢复

三种机制分工不同：

| 机制 | 作用域 | 触发条件 | 典型场景 |
|------|--------|---------|---------|
| `retries:` | rule 级 | 该 rule 执行失败 | 网络闪断（下载参考基因组） |
| `restart-times` | pipeline 级 | 整个 pipeline 失败 | 集群节点故障 |
| `--rerun-incomplete` | 手动 | 用户中断后恢复 | Ctrl+C 终止后继续 |

```yaml
# profile/slurm-apptainer/config.yaml:62 — restart-times
restart-times: 2
```
> 🔗 [profile/slurm-apptainer/config.yaml:62](profile/slurm-apptainer/config.yaml#L62)

### 生命周期钩子

```python
# 成功后发送通知
onsuccess:
    "curl -X POST -d 'pipeline done' https://hooks.slack.com/..."

# 失败后清理
onerror:
    "echo 'Pipeline failed!' | mail -s 'Snakemake Error' admin@lab.org"

# 开始时记录
onstart:
    "date >> pipeline_history.log"
```

### `cache: True` 内容寻址 vs 时间戳缓存

```python
# hisat2.smk:15 — 内容寻址跨批次缓存
cache: True,
```
> 🔗 [hisat2.smk:15](workflow/rules/hisat2.smk#L15)

| 特性 | 时间戳缓存（默认） | 内容寻址（cache: True） |
|------|-------------------|------------------------|
| 判断依据 | mtime | 文件内容哈希 |
| 跨批次复用 | ❌ 不同目录不同时间戳 | ✅ 相同内容即命中 |
| 适用场景 | 常规 rule | 索引构建、参考数据 |

### ⭐ 最佳实践

- **retries 用于瞬态故障**（网络、IO），restart-times 用于节点级故障
- **`--rerun-incomplete` 是中断恢复的首选**：不会重跑已完成的 job
- **cache: True 用于跨批次可复用的产出**：节省索引重建的时间

### ⚠ 常见错误

- **混淆 `retries` 和 `restart-times`**：前者 rule 内重试，后者整个 pipeline 重新提交
- **`--rerun-incomplete` 后不用 `--forcerun`**：两者逻辑不同，不需要同时使用
- **cache: True 时改动 input 内容但保持相同文件名**：缓存基于内容哈希，会正确识别变化

### 💡 面试参考 → 附录 D-18, D-19

---

## 16. 从零搭建 Snakemake 项目

> 跟着分步教程，从空白目录开始搭建一个完整的 Snakemake 工作流。每一步都引用前面章节的知识，最后回顾总结最佳实践。

### 16.1 项目初始化 — 目录结构

从零开始，先创建项目骨架：

```bash
$ mkdir -p my-rnaseq/{workflow/{rules,scripts,envs,data/{raw,ref}},config,profile,containers}
$ cd my-rnaseq
$ touch workflow/Snakefile
$ touch config/config.yaml
```

推荐目录结构（也是本项目的实际布局）：

```
my-rnaseq/
├── config/
│   └── config.yaml              # 唯一配置入口
├── workflow/
│   ├── Snakefile                #   主编排器
│   ├── rules/                   #   规则模块（.smk）
│   ├── scripts/                 #   辅助脚本（Python/R）
│   ├── envs/                    #   Conda 环境定义
│   └── data/                    #   输入数据
│       ├── raw/                 #     原始 FASTQ
│       └── ref/                 #     参考基因组 / GTF / 接头
├── profile/                     # 执行配置
├── containers/                  # 容器构建
├── runs/                        # 运行时产出
│   └── {batch}/
│       ├── results/
│       ├── logs/
│       └── benchmarks/
└── README.md
```

**关键设计**：`workflow/` 子目录包含所有运行时代码，`runs/{batch}/` 按批次隔离产出。这种布局是 Snakemake 官方推荐的标准结构。

### 16.2 第一个 Rule — 从最小到完整模板

先写一个最简单的 rule，验证你的 Snakemake 环境能正常工作：

```python
# workflow/Snakefile — 最小可运行骨架
rule hello:
    output:
        "results/hello.txt",
    shell:
        "echo 'Snakemake works!' > {output}"
```

```bash
$ snakemake -s workflow/Snakefile -j 1
# 输出: results/hello.txt
```

跑通后，逐步添加 components，形成本项目的统一模板：

```python
rule my_tool:
    input:      ...,       # 文件依赖
    output:     ...,       # 产出声明
    params:     ...,       # 非文件参数
    conda:      ...,       # Conda 环境
    container:  ...,       # 容器镜像
    log:        ...,       # 日志输出
    benchmark:  ...,       # 性能记录
    threads:    ...,       # CPU 核数
    shell:      """
        set -euo pipefail
        mkdir -p $(dirname {output.xxx})
        tool_cmd --input {input.xxx} --output {output.xxx} \
            --threads {threads} {params.extra} \
            > {log} 2>&1
        """
```

> 🔗 [hisat2.smk:8-28](workflow/rules/hisat2.smk#L8-L28) — 统一模板的完整示例

### 16.3 引入 config.yaml

不要在 Snakefile 中硬编码路径和参数。第一步重构：把所有可变量提取到 `config.yaml`。

```yaml
# config/config.yaml
paths:
  raw_dir: "workflow/data/raw"
  ref_dir: "workflow/data/ref"
  batch: "20260722_demo"

reference:
  genome_fasta: "workflow/data/ref/genome.fa"
  annotation_gtf: "workflow/data/ref/genes.gtf"
```

```python
# workflow/Snakefile
configfile: "config/config.yaml"

rule some_tool:
    input:
        fasta = config["reference"]["genome_fasta"],
    # ...
```

> 🔗 [Snakefile:70](workflow/Snakefile#L70) — `configfile:` 加载
> 🔗 [config/config.yaml](config/config.yaml) — 完整配置示例

**为后续做准备**：如果你需要批次隔离，在 Snakefile 中加一行 `.format()` 替换即可：

```python
_batch = config["paths"]["batch"]
config["paths"]["outputs"] = {k: v.format(batch=_batch) for k, v in config["paths"]["outputs"].items()}
```

这形成了路径三层抽象的中间层：`config 模板 → format 替换 → os.path.join 拼接`。后续所有代码只引用已格式化的常量。

> 🔗 [Snakefile:73-76](workflow/Snakefile#L73-L76) — 批次占位符替换

### 16.4 添加 Conda 环境

为你的工具创建隔离的 Conda 环境。记住通道顺序规范：

```yaml
# workflow/envs/my_tool.yaml
name: my_tool_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc=0.12.1
```

```python
rule fastqc_raw:
    conda: "../envs/my_tool.yaml"    # 相对于 .smk 文件
    # ...
```

```bash
$ snakemake -s workflow/Snakefile --software-deployment-method conda -j 4
```

环境文件与 rule 放在同一仓库，相对路径引用。`.smk` 在 `workflow/rules/`，env yaml 在 `workflow/envs/`，所以路径为 `../envs/xxx.yaml`。

> 🔗 [envs/hisat2.yaml](workflow/envs/hisat2.yaml) — 环境定义示例
> 🔗 [hisat2.smk:18](workflow/rules/hisat2.smk#L18) — conda 路径写法

### 16.5 模块化拆分 — Snakefile → rules/*.smk

当一个 Snakefile 超过 100 行，就该拆分了。使用 `include:` 将规则分布到多个 `.smk` 模块中：

```python
# workflow/Snakefile — 主编排器（只保留 import、config、常量、rule all、include）
import os
import sys
from common import get_samples, ...

configfile: "config/config.yaml"
SAMPLES = get_samples(config)

OUT_MY_STEP = config["paths"]["outputs"]["my_step"]

rule all:
    input:
        ...

include: "rules/my_tool.smk"
include: "rules/another_tool.smk"
# ... 按 DAG 层级排列
```

```python
# workflow/rules/my_tool.smk — 独立模块
rule my_tool:
    input:   ...
    output:  ...
    conda:   "../envs/my_tool.yaml"
    container: CONTAINER_URI
    log:     os.path.join(LOG_ROOT, "my_tool/{sample}.log")
    benchmark: os.path.join(BENCH_ROOT, "my_tool/{sample}.txt")
    shell:   "set -euo pipefail\nmkdir -p {OUT_MY_STEP}\n..."
```

**关键**：`include:` 是文本级插入。Snakefile 中定义的 `SAMPLES`、`OUT_*`、`LOG_ROOT` 等常量在 `.smk` 中直接可用。

> 🔗 [Snakefile:126-147](workflow/Snakefile#L126-L147) — 7 条 include 按 DAG 层级排列

### 16.6 PE/SE 适配 — input 函数 + rule 级切换

如果你的工具需要同时支持双端和单端数据，在 `.smk` 的 rule 级别用 `if/else` 切换：

```python
# workflow/rules/my_tool.smk
if config["read_pattern"]["mode"] == "single":
    rule my_tool:
        input:
            lambda wildcards: get_se_input(wildcards, config),
        output:
            result=os.path.join(OUT_DIR, "{sample}.result.txt"),
        shell:
            """
            tool --single-read {input.read} --output {output.result}
            """
else:
    rule my_tool:
        input:
            lambda wildcards: get_pe_input(wildcards, config),
        output:
            result=os.path.join(OUT_DIR, "{sample}.result.txt"),
        shell:
            """
            tool --r1 {input.r1} --r2 {input.r2} --output {output.result}
            """
```

input 函数放在 `common.py` 中，通过 λ 闭包传递 config：

```python
# workflow/scripts/common.py
def get_read_inputs(wildcards, config):
    mode = config["read_pattern"]["mode"]
    raw = config["paths"]["raw_dir"]
    if mode == "single":
        return {"read": f"{raw}/{wildcards.sample}.fastq.gz"}
    return {"r1": f"{raw}/{wildcards.sample}_R1.fastq.gz",
            "r2": f"{raw}/{wildcards.sample}_R2.fastq.gz"}
```

切换只需改 config 一个值：`read_pattern.mode: "single"`。无需动 Snakefile 或任何 rule。

> 🔗 [fastqc.smk:10-53](workflow/rules/fastqc.smk#L10-L53) — PE/SE 在 rule 级切换
> 🔗 [common.py:88-95](workflow/scripts/common.py#L88-L95) — 模式感知的 input 函数

### 16.7 容器 + Profile 配置

为每种运行场景创建 Profile，避免每次手动指定参数：

```yaml
# profile/docker/config.yaml
snakefile: "workflow/Snakefile"
software-deployment-method:
  - apptainer
cores: 16
jobs: 4
latency-wait: 30
printshellcmds: true
```

```bash
$ snakemake --profile profile/docker
```

三层 Profile 覆盖全部场景：

| Profile | 场景 | 关键配置 |
|---------|------|---------|
| `profile/docker/` | 本地开发 | Docker 容器，低并行度 |
| `profile/apptainer/` | 单机生产 | Apptainer 容器，满核并行 |
| `profile/slurm-apptainer/` | HPC 集群 | Slurm + Apptainer，资源分级 |

构建容器镜像（推荐统一镜像 — 所有工具装在一个镜像里）：

```bash
$ bash containers/build.sh all     # Docker → SIF
$ bash containers/build.sh test    # 验证工具版本
```

> 🔗 [profile/docker/config.yaml](profile/docker/config.yaml) — Docker profile
> 🔗 [profile/slurm-apptainer/config.yaml](profile/slurm-apptainer/config.yaml) — Slurm profile

### 16.8 测试工作流 — dry-run / stub / E2E

新项目搭好后，分三层验证：

**第一层：语法 + DAG 检查**

```bash
$ snakemake -s workflow/Snakefile --lint       # 语法检查
$ snakemake -s workflow/Snakefile -np          # DAG 预览
$ snakemake -s workflow/Snakefile --dag | dot -Tpdf > dag.pdf   # 可视化
```

**第二层：单样本 stub 测试**

在 rule 中用 `run:` 替代 `shell:` 做 stub（模拟产出，不调真实工具）：

```python
rule my_tool:
    output:
        "results/{sample}.txt",
    run:
        shell("touch {output}")   # stub: 只创建空文件，验证 DAG
```

```bash
$ snakemake -s workflow/Snakefile -j 1 --forcerun my_tool
```

**第三层：小数据 E2E 测试**

用最小测试数据集（如 E. coli 的子集）跑完整流程，验证产出内容：

```bash
$ snakemake -s workflow/Snakefile --software-deployment-method conda -j 4
# 检查产出文件
$ head results/06_deseq2/test_vs_ctrl_significant.csv
$ ls results/07_clusterprofiler/
```

> 🔗 [03. 开发工作流与调试](#03-开发工作流与调试) — 所有调试命令详解

### 16.9 回顾与最佳实践

搭建完成。回顾每一步，总结这个流程为什么这样设计：

#### 路径三层抽象

```
config 模板（{batch} 占位符） → format() 替换 → os.path.join() 拼接
```

改 batch → 所有产出路径自动跟随；改 outputs 子目录名 → 所有代码自动跟随。不硬编码路径，换来的是多批次切换时零改动。

#### 每层设计选择

| 决策 | 选择 | 理由 |
|------|------|------|
| 配置管理 | 单一 config.yaml | 修改参数只改一处 |
| 模块拆分 | include: rules/*.smk | 单文件短小（<100行），模块独立 |
| PE/SE 切换 | rule 级 if/else | DAG 构建时知道完整 input/output |
| 路径拼接 | os.path.join() | 跨平台，不乱写 `/` 或 `\` |
| 中间文件 | temp() | 自动清理，节省磁盘 |
| 关键产出 | protected() | 防止误覆盖（如基因组索引） |
| Shell 安全 | set -euo pipefail | 管道错误不会被静默吞噬 |
| Conda channel | conda-forge→bioconda→defaults | Bioconda 官方推荐顺序 |
| 日志隔离 | 每个 rule 独立 log | 并行运行时日志不交叉 |
| 性能追踪 | benchmark 附加到每个 rule | 优化资源分配的量化依据 |

#### 何时用什么复用方式

| 方式 | 复杂度 | 适用场景 |
|------|--------|---------|
| `wrapper:` | 低 | 标准工具（FastQC、STAR、samtools），一行引用 |
| 自写 rule | 中 | 项目特定的工具调用和参数组合 |
| `module` | 高 | 跨项目复用，Snakemake 8+ 模块化 |

```python
# wrapper 示例 — 一行引用官方预构建规则
rule fastqc:
    wrapper:
        "0.118.2/bio/fastqc"
```

#### 可扩展方向

- **`--report`** — 自动生成 HTML 分析报告，包含运行时间、软件版本、结果摘要
- **`snakemake --cache`** — 远程缓存（S3/GCS），团队共享计算结果
- **`storage:`** — Snakemake 8+ 原生远程存储：`storage: "s3://my-bucket/{sample}.bam"`
- **GitHub Actions CI** — 每次 PR 跑 `--lint` + `-np` dry-run 检查
- **Snakemake module** — 跨项目复用，定义一次 `module my_module:`，多项目 `use rule * from my_module`

#### 核心原则

- ✅ **唯一配置入口**（config.yaml）
- ✅ **批次隔离**（`{batch}` 占位符）
- ✅ **启动前校验**（`validate_config()` + `validate_input_files()`）
- ✅ **每个 rule 独立 log + benchmark**
- ✅ **PE/SE 模式在 rule 级切换**（非 shell 内条件）
- ✅ **Conda channel 顺序固定**（conda-forge → bioconda → defaults）
- ✅ **Shell `set -euo pipefail`**
- ✅ **`temp()` 清理中间文件 + `protected()` 保护索引**
- ✅ **路径三层抽象**（config 模板 → format → os.path.join）
- ✅ **容器化优先**（开发用 conda，生产/CI 用容器保持一致）

#### 一句话总结

> 从项目模板开始，从 `rule all` 反向推导，边写边 `--dag` 可视化验证，跑通小数据再做全量。

### 💡 面试参考 → 附录 D-20

---

# 第四部分：附录

## 附录 A：工具版本速查表

| 工具 | 版本 | Conda 环境 |
|------|------|-----------|
| FastQC | 0.12.1 | envs/fastqc.yaml |
| Trimmomatic | 0.39 | envs/trimmomatic.yaml |
| HISAT2 | 2.2.1 | envs/hisat2.yaml |
| Samtools | 1.18 | envs/hisat2.yaml |
| featureCounts | 2.0.6 | envs/featurecounts.yaml |
| MultiQC | 1.21 | envs/multiqc.yaml |
| R | 4.3.2 | envs/deseq2.yaml |
| DESeq2 | 1.42.0 | envs/deseq2.yaml |
| rtracklayer | 1.62.0 | envs/deseq2.yaml |
| clusterProfiler | 4.10.0 | envs/clusterprofiler.yaml |
| enrichplot | 1.22.0 | envs/clusterprofiler.yaml |

## 附录 B：项目文件导航

| 文件 | 角色 |
|------|------|
| [config/config.yaml](config/config.yaml) | 唯一配置入口 |
| [workflow/Snakefile](workflow/Snakefile) | 主编排器 |
| [workflow/scripts/common.py](workflow/scripts/common.py) | Python 辅助函数（16 个） |
| [workflow/rules/hisat2.smk](workflow/rules/hisat2.smk) | HISAT2 索引 + 比对 |
| [workflow/rules/fastqc.smk](workflow/rules/fastqc.smk) | FastQC 质控（raw + trimmed） |
| [workflow/rules/trimmomatic.smk](workflow/rules/trimmomatic.smk) | Trimmomatic 去接头 + 质量剪切 |
| [workflow/rules/featurecounts.smk](workflow/rules/featurecounts.smk) | featureCounts 基因表达定量 |
| [workflow/rules/multiqc.smk](workflow/rules/multiqc.smk) | MultiQC 汇总质控报告 |
| [workflow/rules/deseq2.smk](workflow/rules/deseq2.smk) | DESeq2 差异表达分析 |
| [workflow/rules/clusterprofiler.smk](workflow/rules/clusterprofiler.smk) | GO + KEGG 功能富集分析 |
| [workflow/scripts/deseq2.R](workflow/scripts/deseq2.R) | DESeq2 R 分析脚本 |
| [workflow/scripts/clusterprofiler.R](workflow/scripts/clusterprofiler.R) | clusterProfiler R 富集脚本 |
| [workflow/scripts/gene2symbol.R](workflow/scripts/gene2symbol.R) | 基因 ID → Symbol 注释 |
| [workflow/scripts/common.R](workflow/scripts/common.R) | R 公共函数（parse_arg） |
| [profile/docker/](profile/docker/) | Docker 本地执行配置 |
| [profile/apptainer/](profile/apptainer/) | Apptainer 单机执行配置 |
| [profile/slurm-apptainer/](profile/slurm-apptainer/) | Slurm 集群执行配置 |
| [containers/](containers/) | 容器构建（Dockerfile + Apptainer + build.sh） |

## 附录 C：常用命令速查表

| 场景 | 命令 |
|------|------|
| Dry-run 预览（含 shell） | `snakemake -s workflow/Snakefile -np` |
| 打印实际 shell 命令 | `snakemake -s workflow/Snakefile -p` |
| 文件状态表 | `snakemake -s workflow/Snakefile --summary` |
| 强制重跑后状态表 | `snakemake -s workflow/Snakefile --summary --forcerun <rule>` |
| DAG 可视化 | `snakemake -s workflow/Snakefile --dag \| dot -Tpdf > dag.pdf` |
| 强制重跑某 rule | `snakemake -s workflow/Snakefile --forcerun <rule>` |
| 中断恢复 | `snakemake -s workflow/Snakefile --rerun-incomplete` |
| 生成 HTML 报告 | `snakemake -s workflow/Snakefile --report report.html` |
| Conda 模式运行 | `snakemake -s workflow/Snakefile --software-deployment-method conda -j 8` |
| Docker 模式运行 | `snakemake -s workflow/Snakefile --profile profile/docker` |
| Apptainer 模式运行 | `snakemake -s workflow/Snakefile --profile profile/apptainer` |
| Slurm 集群模式 | `snakemake -s workflow/Snakefile --profile profile/slurm-apptainer` |
| 语法检查 | `snakemake -s workflow/Snakefile --lint` |
| 标记为已完成 | `snakemake -s workflow/Snakefile --touch <rule>` |
| 指定调度策略 | `snakemake -s workflow/Snakefile --scheduler ilp\|greedy` |
| 单独构建某目标 | `snakemake -s workflow/Snakefile <target_file_path>` |

## 附录 D：知识自测 20 题（面试常见）

每题为 "Q → A → ❌ 常见误区" 格式，标注关联章节号。

### D-1（01 DAG）逆向推导

**Q**：`snakemake -s workflow/Snakefile` 不指定目标，会执行哪些 rule？

**A**：Snakefile 中**文本顺序上第一个出现的 `rule` 定义**作为默认目标（`include:` 是文本级插入，展开后按整体文本顺序算）。本项目 `rule all`（第 114 行）在所有 `include:` 之前，是第一个 rule，Snakemake 从它的 `input` 列表逆向推导，执行所有依赖链上的 rule。

**❌ 常见误区**：① 以为不指定目标会执行所有 rule（实际只执行默认目标的依赖链）；② 以为 include 的 rule 不算"第一个"（include 是文本级插入，无边界）。例如本项目若删掉 `rule all`，默认目标变成 `include` 展开后的 `rule hisat2_index`，只会构建索引，比对、定量、差异分析全部跳过。

### D-2（01 DAG）并行度判断

**Q**：DAG 中有 20 个 job，`-j 4` 时多少个同时运行？

**A**：取决于 DAG 层级。同层无依赖的 job 最多 4 个同时运行（受 `-j` 限制），不同层的 job 按依赖关系串行。如果 Layer 1 有 6 个无依赖 job，最多 4 个并行，剩余 2 个等槽位释放。

**❌ 常见误区**：以为 `-j 4` 就是 4 个 job 从头跑到尾。实际是 DAG 每层最多 4 个并行。

### D-3（02 Rule）shell 健壮性

**Q**：为什么 Shell 块首行要写 `set -euo pipefail`？

**A**：
- `-e`：任何命令失败（退出码非 0）立即退出
- `-u`：引用未定义变量时报错
- `-o pipefail`：管道中任一命令失败，整个管道视为失败

在 `hisat2 | samtools sort` 这种管道中，没有 `pipefail` 的话 hisat2 失败时 shell 只检查管道最后一个命令（samtools）的退出码，静默吞噬错误。

**❌ 常见误区**：以为 `-e` 足够——`-e` 不处理管道中间命令的失败。

### D-4（02 Rule）input 函数 vs 静态列表

**Q**：什么时候用 input 函数，什么时候用静态列表？

**A**：
- 静态列表：所有样本的 input 来源和格式一致（如 `config["reference"]["genome_fasta"]`）
- input 函数：input 路径依赖 `{sample}` wildcard 动态生成（如 `get_read_input_list(wildcards, config)`）

**❌ 常见误区**：在 input 函数外面使用 `wildcards` 变量——wildcards 只在函数闭包内可用。

### D-5（04 Wildcards）名称一致性

**Q**：为什么 output 和 input 中 wildcard 名称必须一致？

**A**：Snakemake 通过匹配 output 和 input 中同名 wildcard 来推导值。如 output 中 `{sample}` 匹配到 `WT_1`，则 input 中 `{sample}` 也替换为 `WT_1`。名称不一致会导致 Snakemake 无法推导 wildcard 值。

**❌ 常见误区**：以为 wildcard 是全局变量——实际上每个 rule 独立推导。

### D-6（04 Wildcards）glob_wildcards vs 手动列表

**Q**：`glob_wildcards()` 相比手动列出样本的优劣？

**A**：优势是自动发现（无需维护样本表），劣势是依赖文件系统状态（文件被误删会导致样本丢失）。推荐：开发探索用 `glob_wildcards()`，生产用 `samples_file`（CSV/TSV）明确管理。

**❌ 常见误区**：`glob_wildcards()` 读取文件系统，不会读取 config 或样本表。

### D-7（05 expand）output 不能用 expand

**Q**：为什么 `output:` 不能使用 `expand()`？

**A**：`expand()` 生成完整列表（笛卡尔积），而 `output:` 应用 `{wildcard}` 声明本 rule **一次运行**的产出。Snakemake 需要从 output 路径中提取 wildcard 值来驱动 job 调度。如果用 expand，Snakemake 无法知道哪些文件由哪个 wildcard 值产生。

**❌ 常见误区**：以为 expand 在 output 和 input 中用法相同。

### D-8（05 expand）多维 expand 笛卡尔积陷阱

**Q**：`expand("data/{sample}_{read}.fq", sample=["A","B"], read=["R1","R2"])` 生成什么？

**A**：4 个文件：`A_R1.fq, A_R2.fq, B_R1.fq, B_R2.fq`（笛卡尔积）。如果想一一对应（`A_R1.fq, B_R2.fq`），需要 `expand("data/{sample}_{read}.fq", zip, sample=["A","B"], read=["R1","R2"])`。

**❌ 常见误区**：不检查 expand 生成的路径数量，导致 input 收集到不存在的文件或遗漏存在的文件。

### D-9（06 Config）Python 层 vs shell 层引用

**Q**：`config["key"]` 在 Python 层和 shell 层的引用方式有什么区别？

**A**：
- Python 层（input/output/params）：`config["reference"]["genome_fasta"]`
- Shell 层：`{config[reference][genome_fasta]}`（不用引号，用方括号）
- Params 中转后在 shell 层：`{params.prefix}`

**❌ 常见误区**：在 shell 中写 `{config["key"]}` 带引号，这会导致 shell 解析错误。

### D-10（06 Config）params vs resources

**Q**：params 和 resources 的分工是什么？

**A**：
- `params:` — 工具参数（传给程序的命令行参数或配置值）
- `resources:` — 资源限制（mem_mb / runtime / cpus_per_task / tmpdir），传给调度器
- `threads:` — CPU 核数，是 resources 的简化版

**❌ 常见误区**：把内存限制放在 params 里（不会限制实际资源使用），或把命令行参数放在 resources 里（不会传给程序）。

### D-11（07 Conda）conda 与 container 为什么互斥

**Q**：为什么 `--software-deployment-method conda` 与 container 互斥？

**A**：两者是替代关系。conda 模式在运行时创建隔离环境安装依赖；container 模式使用预构建镜像（已包含全部工具）。同时启用会导致：容器内的基础环境 + conda 环境叠加，路径冲突、版本混乱。开发用 conda（快速迭代），生产用 container（一致性部署）。

**❌ 常见误区**：以为 conda + container 可以互补——实际上容器镜像已经包含所有工具，不需要 conda。

### D-12（11 temp/protected）temp() 删除时机

**Q**：`temp()` 标记的文件什么时候被删除？

**A**：当所有以该文件为 input 的下游 rule 都成功完成后，Snakemake 自动删除。如果某个下游 rule 还在等待执行或失败重试中，temp 文件保留。

**❌ 常见误区**：以为 temp 文件在当前 rule 完成后立即删除——需要等下游全部完成。

### D-13（11 temp/protected）protected() 如何覆盖

**Q**：`protected()` 保护的文件如何强制重跑？

**A**：`snakemake --forceall --protected`。需要同时指定两个 flag：`--forceall`（忽略时间戳）、`--protected`（允许覆盖写保护文件）。单独使用 `--forceall` 会被 `protected()` 阻挡。

**❌ 常见误区**：以为 `--forcerun` 或 `-f` 就能覆盖 protected 文件。

### D-14（12 ruleorder）歧义消除

**Q**：两个 rule 都能产生同一个 output 文件时，Snakemake 的行为是什么？

**A**：Snakemake 抛出 `AmbiguousRuleException`，不会自动选择。需要用 `ruleorder: A > B` 消除歧义，或者重新设计 output 路径使两个 rule 产出不同的文件。

**❌ 常见误区**：以为 Snakemake 会按某种规则（如 rule 定义顺序）自动选择——实际上直接报错。

### D-15（13 checkpoint）核心区别

**Q**：checkpoint 和普通 rule 的核心区别是什么？

**A**：普通 rule 的 output 在 DAG 构建阶段完全确定（静态），checkpoint 的 output 在运行时才知道。checkpoint 运行后可能产生新的 wildcard 值，下游规则必须通过 `checkpoints.rule_name.get()` 获取运行结果后才能构建 DAG。

**❌ 常见误区**：用 checkpoint 解决可以用 input 函数解决的问题——增加了不必要的复杂度。

### D-16（14 HPC）jobs vs cores

**Q**：`jobs` 和 `cores` 的区别是什么？

**A**：
- `cores`：单机可用总核数。Snakemake 将任务调度到本地 CPU。每个 job 的 `threads` 从这个池中分配。
- `jobs`：集群模式允许同时提交到 Slurm 的最大作业数。每个作业独立申请 `cpus_per_task` 核。

本地模式用 `cores`，集群模式（`executor: slurm`）用 `jobs`。

**❌ 常见误区**：在 Slurm 配置中设置 `cores` 期望它能限制集群核数——集群核数由 Slurm 分配，`jobs` 只限制同时提交数。

### D-17（14 HPC）latency-wait + localrules

**Q**：为什么 NFS 环境需要 `latency-wait`？`localrules:` 的作用是什么？

**A**：
- `latency-wait`：NFS 文件系统有写入延迟。Job A 完成并写出 output 文件后，Job B 所在的节点可能暂时看不到该文件。`latency-wait` 让 Snakemake 等待指定秒数再检查文件存在性。
- `localrules:`：标记在登录节点（而非计算节点）运行的规则。`rule all`、目录创建等秒级任务不应提交到 Slurm 排队。

**❌ 常见误区**：以为 `latency-wait` 只在集群需要——任何网络文件系统（NFS/Lustre/GPFS）都可能需要。

### D-18（15 生命周期）失败处理三分工

**Q**：`retries`、`restart-times`、`--rerun-incomplete` 三者的分工？

**A**：
- `retries: 3`（rule 级）：同一 rule 执行失败后立即重试 3 次，适合网络闪断
- `restart-times: 2`（pipeline 级）：整个 pipeline 失败后重试 2 次，适合集群节点故障
- `--rerun-incomplete`（手动）：用户 Ctrl+C 中断后，下次运行时加上此参数恢复

**❌ 常见误区**：用 `restart-times` 替代 `retries`——restart-times 重跑整个 pipeline 代价大，retries 只重跑失败的单个 rule。

### D-19（15 生命周期）cache: True 内容寻址

**Q**：`cache: True` 内容寻址缓存和时间戳缓存的核心区别？

**A**：
- 时间戳（默认）：比较 input 的 mtime 是否新于 output → 决定是否重跑
- 内容寻址（`cache: True`）：基于文件内容哈希判断 → 相同内容的 input 产生相同 output → 跨批次复用

`cache: True` 适合索引构建：同一个基因组 fasta 文件，不同批次间即使路径不同，只要内容相同就能命中缓存。

**❌ 常见误区**：以为 `cache: True` 基于文件路径缓存——实际基于内容的哈希值。

### D-20（16 从零搭建）路径多层抽象

**Q**：为什么要做 config 模板 → format → os.path.join 三层路径抽象，而不是直接硬编码路径？

**A**：
1. **config 模板层**（`{batch}` 占位符）：批次数值集中管理，改一处全量切换
2. **format 替换层**（Snakefile 中 `.format(batch=...)`）：运行时动态填充，Python 代码中统一处理
3. **os.path.join 层**（.smk 规则中）：跨平台兼容，不需手动处理 `/` vs `\`

三层设计实现：改 batch → 所有产出路径自动跟随；改 outputs 子目录名 → 所有代码自动跟随。

**❌ 常见误区**：以为加一层抽象是过度设计——在一个多批次项目中，硬编码路径的维护成本会指数增长。

---

## 附录 E：DAG 拓扑图（本项目）

```
Layer 0:  hisat2_index ─────────────────────────────────────┐
                                                             │
Layer 1:  fastqc_raw ─────────────────────────────────────┐ │
                                                           │ │
Layer 2:  trimmomatic ──> fastqc_trimmed ──┐              │ │
              │  │                           │              │ │
              │  └──> hisat2_align ──────────┤──> multiqc   │ │
              │        (Layer 3b)            │   (Layer 5)  │ │
              │                              │              │ │
Layer 4:  featureCounts <───────────────────┘              │ │
              │                                              │ │
Layer 6:  DESeq2                                            │ │
              │                                              │ │
Layer 7:  clusterProfiler                                   │ │
```

| Layer | Rule | 并行度 | 说明 |
|-------|------|--------|------|
| 0 | hisat2_index | 1 | 一次性建索引，无 wildcard |
| 1 | fastqc_raw | N | 每个 `{sample}` 一次 |
| 2 | trimmomatic | N | 每个 `{sample}` 一次，依赖 Layer 1 |
| 3 | fastqc_trimmed | N | 每个 `{sample}` 一次，依赖 Layer 2 |
| 3 | hisat2_align | N | 每个 `{sample}` 一次，依赖 Layer 2 + Layer 0 |
| 4 | featureCounts | 1 | 汇总所有 BAM |
| 5 | MultiQC | 1 | 汇总所有 QC 报告 |
| 6 | DESeq2 | 1 | 依赖 Layer 4 |
| 7 | clusterProfiler | 1 | 依赖 Layer 6 |

> 🔗 [README_CN.md#流程-dag](README_CN.md#流程-dag) — 流程 DAG 图解
> 🔗 [Snakefile:1-34](workflow/Snakefile#L1-L34) — DAG 概览注释 + 模块映射
