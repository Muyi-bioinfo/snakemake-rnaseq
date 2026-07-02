# ============================================================
# clusterProfiler 模块 — GO + KEGG 功能富集分析（下游）
# 依赖: deseq2 (显著差异基因列表)
#
# 输入: 06_deseq2/{contrast}_significant.csv
# 输出: 07_clusterprofiler/{contrast}_GO/KEGG_enrichment.csv + *dotplot.pdf
# ============================================================

rule clusterprofiler:
    input:
        # 收集所有对比的显著基因 CSV
        sig_files=lambda wildcards: [
            os.path.join(OUT_DESEQ2, f"{c['name']}_significant.csv")
            for c in config["deseq2"]["contrasts"]
        ],
    output:
        *get_clusterprofiler_outputs(config),
    conda:
        "../envs/clusterprofiler.yaml"
    container:
        CONTAINER_URI
    log:
        os.path.join(LOG_ROOT, "clusterprofiler.log"),
    benchmark:
        os.path.join(BENCH_ROOT, "clusterprofiler.txt"),
    params:
        contrast_names=lambda wildcards: ",".join(
            c["name"] for c in config["deseq2"]["contrasts"]
        ),
        sig_files_str=lambda wildcards: ",".join(
            os.path.join(OUT_DESEQ2, f"{c['name']}_significant.csv")
            for c in config["deseq2"]["contrasts"]
        ),
        org_db=config["clusterprofiler"]["org_db"],
        kegg_org=config["clusterprofiler"]["kegg_organism"],
        from_type=config["clusterprofiler"]["from_type"],
        pval=config["clusterprofiler"]["pvalue_cutoff"],
        qval=config["clusterprofiler"]["qvalue_cutoff"],
        gene_id_col=config["clusterprofiler"]["gene_id_col"],
        show_cat=config["clusterprofiler"]["show_category"],
        outdir=config["paths"]["outputs"]["clusterprofiler"],
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        Rscript workflow/scripts/clusterprofiler.R \
            --sig_files {params.sig_files_str} \
            --contrast_names {params.contrast_names} \
            --org_db {params.org_db} \
            --kegg_org {params.kegg_org} \
            --from_type {params.from_type} \
            --pval {params.pval} \
            --qval {params.qval} \
            --gene_id_col {params.gene_id_col} \
            --show_cat {params.show_cat} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """
