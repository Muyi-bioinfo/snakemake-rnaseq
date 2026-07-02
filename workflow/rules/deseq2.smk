# ============================================================
# DESeq2 模块 — 差异表达分析（下游）
# 依赖: featurecounts (计数矩阵)
#
# 输入: 04_featurecounts/featurecounts.txt
# 输出: 06_deseq2/*.csv + *.pdf
# ============================================================

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
    benchmark:
        os.path.join(BENCH_ROOT, "deseq2.txt"),
    params:
        groups_json=get_deseq2_groups_json(config),
        contrasts_json=get_deseq2_contrasts_json(config),
        outdir=config["paths"]["outputs"]["deseq2"],
        padj=config["deseq2"]["padj_threshold"],
        log2fc=config["deseq2"]["log2fc_threshold"],
        top_n=config["deseq2"]["top_n_genes"],
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
            --top_n {params.top_n} \
            > {log} 2>&1
        """
