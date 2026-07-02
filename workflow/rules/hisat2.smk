# ============================================================
# HISAT2 模块 — 索引构建 + 序列比对
# 依赖: reference.fasta (索引构建), trimmomatic (比对)
# ============================================================

# ── Layer 0: hisat2_index ──────────────────────────────────

rule hisat2_index:
    input:
        fasta=config["reference"]["genome_fasta"],
    output:
        protected(touch("workflow/data/ref/hisat2_index/genome.1.ht2")),
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


# ── Layer 3b: hisat2_align ─────────────────────────────────
# PE/SE 在 rule 级切换：input key 名称不同，Snakemake shell 模板要求所有引用键都存在

if config["read_pattern"]["mode"] == "single":

    rule hisat2_align:
        input:
            unpack(lambda wildcards: get_hisat2_input(wildcards, config)),
        output:
            bam=os.path.join(OUT_HISAT2, "{sample}.sorted.bam"),
            bai=os.path.join(OUT_HISAT2, "{sample}.sorted.bam.bai"),
        conda: "../envs/hisat2.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "hisat2_align/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "hisat2_align/{sample}.txt")
        threads: config["hisat2"]["threads"]
        params:
            index_prefix=config["reference"]["hisat2_index_prefix"],
            extra=config["hisat2"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][hisat2]}
            hisat2 -x {params.index_prefix} \
                -U {input.read} \
                -p {threads} {params.extra} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
            """

else:

    rule hisat2_align:
        input:
            unpack(lambda wildcards: get_hisat2_input(wildcards, config)),
        output:
            bam=os.path.join(OUT_HISAT2, "{sample}.sorted.bam"),
            bai=os.path.join(OUT_HISAT2, "{sample}.sorted.bam.bai"),
        conda: "../envs/hisat2.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "hisat2_align/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "hisat2_align/{sample}.txt")
        threads: config["hisat2"]["threads"]
        params:
            index_prefix=config["reference"]["hisat2_index_prefix"],
            extra=config["hisat2"]["extra"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][hisat2]}
            hisat2 -x {params.index_prefix} \
                -1 {input.r1} -2 {input.r2} \
                -p {threads} {params.extra} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
            """
