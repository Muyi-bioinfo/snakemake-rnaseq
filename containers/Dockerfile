# ============================================================
# RNA-seq 上游分析流程 — 统合容器镜像
# ============================================================
# 包含全部工具：FastQC, Trimmomatic, HISAT2, Samtools,
#               featureCounts, MultiQC, R + DESeq2 + clusterProfiler
#
# 使用方式:
#   docker build -t rnaseq-pipeline:latest -f containers/Dockerfile .
#
# 验证:
#   docker run --rm rnaseq-pipeline:latest fastqc --version
#   docker run --rm rnaseq-pipeline:latest Rscript -e "library(DESeq2)"
#   docker run --rm rnaseq-pipeline:latest Rscript -e "library(clusterProfiler)"
# ============================================================

FROM continuumio/miniconda3:24.1.2-0

# ---- 系统依赖 ----
RUN apt-get update -qq && \
    apt-get install -y -qq --no-install-recommends \
        procps \
        curl \
        wget \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ---- 配置 conda channels ----
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# ---- 安装全部生物信息学工具 ----
# 分两层安装：先装 R 生态（耗时较长），再装命令行工具
# 这样利于缓存复用

# Layer 1: R + Bioconductor 包
RUN conda install -y -q \
        r-base=4.3.2 \
        bioconductor-deseq2=1.42.0 \
        r-ggplot2=3.5.0 \
        r-pheatmap=1.0.12 \
        r-ggrepel=0.9.5 \
        r-readr=2.1.5 \
        r-dplyr=1.1.4 \
        r-tidyr=1.3.1 \
        r-jsonlite=1.8.8 \
        bioconductor-rtracklayer=1.62.0 \
        bioconductor-clusterprofiler=4.10.0 \
        bioconductor-enrichplot=1.22.0 \
        bioconductor-org.hs.eg.db=3.18.0 \
        bioconductor-org.mm.eg.db=3.18.0 \
        r-cowplot=1.1.3 \
        && \
    conda clean -afy

# Layer 2: 命令行工具
RUN conda install -y -q \
        fastqc=0.12.1 \
        trimmomatic=0.39 \
        hisat2=2.2.1 \
        samtools=1.18 \
        subread=2.0.6 \
        multiqc=1.21 \
        && \
    conda clean -afy

# ---- 环境变量 ----
ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# ---- 创建标准目录结构 ----
RUN mkdir -p /pipeline/workflow/{data/{raw,ref},scripts} && \
    mkdir -p /pipeline/{results,logs,benchmarks} && \
    chmod -R 777 /pipeline

WORKDIR /pipeline

# ---- 元数据 ----
LABEL org.opencontainers.image.title="RNA-seq Pipeline" \
      org.opencontainers.image.description="Snakemake RNA-seq upstream analysis: FastQC → Trimmomatic → HISAT2 → featureCounts → MultiQC → DESeq2 → clusterProfiler" \
      org.opencontainers.image.authors="pipeline-team" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.source="https://github.com/example/rnaseq-pipeline"

# ---- 默认命令 ----
ENTRYPOINT ["/bin/bash"]
CMD ["snakemake", "--help"]
