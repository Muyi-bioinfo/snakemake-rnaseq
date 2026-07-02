# ============================================================
# Trimmomatic 模块 — 去接头 + 质量剪切
# 依赖: fastqc_raw
#
# PE/SE 在 rule 级 + shell 级双重切换
# ============================================================

if config["read_pattern"]["mode"] == "single":

    rule trimmomatic:
        input:
            unpack(lambda wildcards: get_read_inputs(wildcards, config)),
        output:
            read=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}.trimmed.fastq.gz")),
        conda: "../envs/trimmomatic.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "trimmomatic/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "trimmomatic/{sample}.txt")
        threads: config["trimmomatic"]["threads"]
        params:
            adapter=config["trimmomatic"]["adapter_fasta"],
            illuminaclip=config["trimmomatic"]["illuminaclip"],
            leading=config["trimmomatic"]["leading"],
            trailing=config["trimmomatic"]["trailing"],
            slidingwindow=config["trimmomatic"]["slidingwindow"],
            minlen=config["trimmomatic"]["minlen"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][trimmomatic]}
            trimmomatic SE -threads {threads} \
                {input.read} {output} \
                ILLUMINACLIP:{params.adapter}:{params.illuminaclip} \
                LEADING:{params.leading} TRAILING:{params.trailing} \
                SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} \
                > {log} 2>&1
            """

else:

    rule trimmomatic:
        input:
            unpack(lambda wildcards: get_read_inputs(wildcards, config)),
        output:
            r1=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R1.trimmed.fastq.gz")),
            r2=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R2.trimmed.fastq.gz")),
            r1_unpaired=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R1.unpaired.fastq.gz")),
            r2_unpaired=temp(os.path.join(OUT_TRIMMOMATIC, "{sample}_R2.unpaired.fastq.gz")),
        conda: "../envs/trimmomatic.yaml"
        container: CONTAINER_URI
        log:     os.path.join(LOG_ROOT, "trimmomatic/{sample}.log")
        benchmark: os.path.join(BENCH_ROOT, "trimmomatic/{sample}.txt")
        threads: config["trimmomatic"]["threads"]
        params:
            adapter=config["trimmomatic"]["adapter_fasta"],
            illuminaclip=config["trimmomatic"]["illuminaclip"],
            leading=config["trimmomatic"]["leading"],
            trailing=config["trimmomatic"]["trailing"],
            slidingwindow=config["trimmomatic"]["slidingwindow"],
            minlen=config["trimmomatic"]["minlen"],
        shell:
            """
            set -euo pipefail
            mkdir -p {config[paths][outputs][trimmomatic]}
            trimmomatic PE -threads {threads} \
                {input.r1} {input.r2} \
                {output.r1} {output.r1_unpaired} \
                {output.r2} {output.r2_unpaired} \
                ILLUMINACLIP:{params.adapter}:{params.illuminaclip} \
                LEADING:{params.leading} TRAILING:{params.trailing} \
                SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} \
                > {log} 2>&1
            """
