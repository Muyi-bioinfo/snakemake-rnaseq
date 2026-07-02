# ============================================================
# MultiQC 模块 — 汇总质控报告
# 依赖: fastqc_raw, fastqc_trimmed, hisat2_align (日志)
#
# 扫描 results/ + logs/ 目录，输出到 05_multiqc/
# ============================================================

rule multiqc:
    input:
        fastqc_raw=expand(
            [os.path.join(OUT_FASTQC_RAW, "{sample}_fastqc.html"),
             os.path.join(OUT_FASTQC_RAW, "{sample}_fastqc.zip")],
            sample=SAMPLES,
        ) if config["read_pattern"]["mode"] == "single" else expand(
            [os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.html"),
             os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.html"),
             os.path.join(OUT_FASTQC_RAW, "{sample}_R1_fastqc.zip"),
             os.path.join(OUT_FASTQC_RAW, "{sample}_R2_fastqc.zip")],
            sample=SAMPLES,
        ),
        fastqc_trimmed=expand(
            [os.path.join(OUT_FASTQC_TRIMMED, "{sample}.trimmed_fastqc.html"),
             os.path.join(OUT_FASTQC_TRIMMED, "{sample}.trimmed_fastqc.zip")],
            sample=SAMPLES,
        ) if config["read_pattern"]["mode"] == "single" else expand(
            [os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R1.trimmed_fastqc.html"),
             os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R2.trimmed_fastqc.html"),
             os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R1.trimmed_fastqc.zip"),
             os.path.join(OUT_FASTQC_TRIMMED, "{sample}_R2.trimmed_fastqc.zip")],
            sample=SAMPLES,
        ),
        hisat2_logs=expand(f"{LOG_ROOT}/hisat2_align/{{sample}}.log", sample=SAMPLES),
    output:
        os.path.join(OUT_MULTIQC, "multiqc_report.html"),
    conda:
        "../envs/multiqc.yaml"
    container:
        CONTAINER_URI
    log:
        os.path.join(LOG_ROOT, "multiqc.log"),
    benchmark:
        os.path.join(BENCH_ROOT, "multiqc.txt"),
    shell:
        """
        set -euo pipefail
        mkdir -p {config[paths][outputs][multiqc]}

        multiqc $(dirname {config[paths][outputs][multiqc]}) {config[paths][log_root]} \
            --outdir {config[paths][outputs][multiqc]} \
            --force \
            > {log} 2>&1
        """
