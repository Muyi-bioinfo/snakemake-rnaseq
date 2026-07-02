# ============================================================
# featureCounts 模块 — 基因表达定量
# 依赖: hisat2_align (所有样本的 BAM)
#
# 输入: 03_hisat2/{sample}.sorted.bam
# 输出: 04_featurecounts/featurecounts.txt + .summary.txt
# ============================================================

rule featurecounts:
    input:
        bams=[f"{OUT_HISAT2}/{s}.sorted.bam" for s in SAMPLES],
        gtf=config["reference"]["annotation_gtf"],
    output:
        counts=os.path.join(OUT_FEATURECOUNTS, "featurecounts.txt"),
        summary=os.path.join(OUT_FEATURECOUNTS, "featurecounts.summary.txt"),
    conda:
        "../envs/featurecounts.yaml"
    container:
        CONTAINER_URI
    log:
        os.path.join(LOG_ROOT, "featurecounts.log"),
    benchmark:
        os.path.join(BENCH_ROOT, "featurecounts.txt"),
    threads: config["featurecounts"]["threads"]
    params:
        feature_type=config["featurecounts"]["feature_type"],
        attr_type=config["featurecounts"]["attr_type"],
        strandedness=config["featurecounts"]["strandedness"],
        extra=config["featurecounts"]["extra"],
    shell:
        """
        set -euo pipefail
        mkdir -p {config[paths][outputs][featurecounts]}

        featureCounts -T {threads} \
            -t {params.feature_type} \
            -g {params.attr_type} \
            -s {params.strandedness} \
            -a {input.gtf} \
            -o {output.counts} \
            {params.extra} \
            {input.bams} \
            > {log} 2>&1

        mv {config[paths][outputs][featurecounts]}/featurecounts.txt.summary {output.summary}
        """
