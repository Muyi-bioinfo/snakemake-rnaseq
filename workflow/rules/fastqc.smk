# ============================================================
# FastQC 模块 — 原始数据 + 剪切后数据质控
# 依赖: data/raw (raw), trimmomatic (trimmed)
#
# PE/SE 在 rule 级切换：输出文件数量不同，Snakemake 9 不支持 output 函数
# ============================================================

# ── Layer 1: fastqc_raw ────────────────────────────────────

if config["read_pattern"]["mode"] == "single":
    rule fastqc_raw:
        input:
            lambda wildcards: get_read_input_list(wildcards, config),
        output:
            html=os.path.join(OUT_FASTQC_RAW, "{sample}_fastqc.html"),
            zip=os.path.join(OUT_FASTQC_RAW,  "{sample}_fastqc.zip"),
        conda: "../envs/fastqc.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "fastqc_raw/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "fastqc_raw/{sample}.txt")
        threads: config["fastqc"]["threads"]
        params:
            extra=config["fastqc"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][fastqc_raw]}
            fastqc {input} --outdir {config[paths][outputs][fastqc_raw]} \
                --threads {threads} {params.extra} > {log} 2>&1
            """
else:
    rule fastqc_raw:
        input:
            lambda wildcards: get_read_input_list(wildcards, config),
        output:
            html_r1=os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.html"),
            html_r2=os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.html"),
            zip_r1=os.path.join(OUT_FASTQC_RAW,  "{sample}_R1_fastqc.zip"),
            zip_r2=os.path.join(OUT_FASTQC_RAW,  "{sample}_R2_fastqc.zip"),
        conda: "../envs/fastqc.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "fastqc_raw/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "fastqc_raw/{sample}.txt")
        threads: config["fastqc"]["threads"]
        params:
            extra=config["fastqc"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][fastqc_raw]}
            fastqc {input} --outdir {config[paths][outputs][fastqc_raw]} \
                --threads {threads} {params.extra} > {log} 2>&1
            """


# ── Layer 3a: fastqc_trimmed ───────────────────────────────

if config["read_pattern"]["mode"] == "single":
    rule fastqc_trimmed:
        input:
            lambda wildcards: get_trimmed_input_list(wildcards, config),
        output:
            html=os.path.join(OUT_FASTQC_TRIMMED, "{sample}.trimmed_fastqc.html"),
            zip=os.path.join(OUT_FASTQC_TRIMMED,  "{sample}.trimmed_fastqc.zip"),
        conda: "../envs/fastqc.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "fastqc_trimmed/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "fastqc_trimmed/{sample}.txt")
        threads: config["fastqc"]["threads"]
        params:
            extra=config["fastqc"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][fastqc_trimmed]}
            fastqc {input} --outdir {config[paths][outputs][fastqc_trimmed]} \
                --threads {threads} {params.extra} > {log} 2>&1
            """
else:
    rule fastqc_trimmed:
        input:
            lambda wildcards: get_trimmed_input_list(wildcards, config),
        output:
            html_r1=os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R1.trimmed_fastqc.html"),
            html_r2=os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R2.trimmed_fastqc.html"),
            zip_r1=os.path.join(OUT_FASTQC_TRIMMED,  "{sample}_R1.trimmed_fastqc.zip"),
            zip_r2=os.path.join(OUT_FASTQC_TRIMMED,  "{sample}_R2.trimmed_fastqc.zip"),
        conda: "../envs/fastqc.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "fastqc_trimmed/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "fastqc_trimmed/{sample}.txt")
        threads: config["fastqc"]["threads"]
        params:
            extra=config["fastqc"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][fastqc_trimmed]}
            fastqc {input} --outdir {config[paths][outputs][fastqc_trimmed]} \
                --threads {threads} {params.extra} > {log} 2>&1
            """
