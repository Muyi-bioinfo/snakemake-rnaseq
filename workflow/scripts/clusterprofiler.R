#!/usr/bin/env Rscript
# ============================================================
# clusterProfiler — 基因功能富集分析（GO + KEGG）
# 输入：DESeq2 显著差异基因 CSV
# 输出：GO/KEGG 富集表格 + 气泡图
#
# 调用方式：
#   Rscript clusterprofiler.R \
#       --sig_files <csv1,csv2,...> \
#       --contrast_names <name1,name2,...> \
#       --org_db org.Hs.eg.db \
#       --kegg_org hsa \
#       --from_type ENSEMBL \
#       --pval 0.05 --qval 0.2 \
#       --gene_id_col gene_id \
#       --show_cat 15 \
#       --outdir results/07_clusterprofiler
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(cowplot)
})

# ---- 参数解析 ----
source("workflow/scripts/common.R")

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("Usage: Rscript clusterprofiler.R --sig_files <csv,...> \\
       --contrast_names <name,...> --org_db <pkg> --kegg_org <code> \\
       --pval <float> --qval <float> --gene_id_col <col> \\
       --show_cat <int> --outdir <dir>\n")
  quit(status = 1)
}

if (length(args) == 0) usage()

sig_files_str  <- parse_arg(args, "--sig_files")
contrast_names_str <- parse_arg(args, "--contrast_names")
org_db_str     <- parse_arg(args, "--org_db")
kegg_org       <- parse_arg(args, "--kegg_org")
from_type      <- parse_arg(args, "--from_type")
pval_cut       <- as.numeric(parse_arg(args, "--pval"))
qval_cut       <- as.numeric(parse_arg(args, "--qval"))
gene_id_col    <- parse_arg(args, "--gene_id_col")
show_cat       <- as.integer(parse_arg(args, "--show_cat"))
outdir         <- parse_arg(args, "--outdir")

# 默认值
if (is.null(pval_cut))   pval_cut   <- 0.05
if (is.null(qval_cut))   qval_cut   <- 0.2
if (is.null(gene_id_col)) gene_id_col <- "gene_id"
if (is.null(show_cat))   show_cat   <- 15
if (is.null(outdir))     outdir     <- "results/07_clusterprofiler"

# ---- 加载 OrgDb ----
if (is.null(org_db_str)) {
  stop("必须通过 --org_db 指定物种注释包（如 org.Hs.eg.db）")
}
suppressPackageStartupMessages({
  library(org_db_str, character.only = TRUE)
})

# 用 get() 获取 OrgDb 对象（clusterProfiler 需要）
org_db <- get(org_db_str)

# ---- 分割逗号分隔的参数 ----
sig_files      <- strsplit(sig_files_str, ",")[[1]]
contrast_names <- strsplit(contrast_names_str, ",")[[1]]

if (length(sig_files) != length(contrast_names)) {
  stop(sprintf("sig_files 数量(%d) 与 contrast_names 数量(%d) 不匹配",
               length(sig_files), length(contrast_names)))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[clusterProfiler] OrgDb: %s | KEGG: %s | pval<=%.3f qval<=%.3f | 展示 top %d",
                org_db_str, kegg_org, pval_cut, qval_cut, show_cat))
message(sprintf("[clusterProfiler] 处理 %d 个对比", length(contrast_names)))

# ============================================================
# 辅助函数
# ============================================================

enrich_and_plot <- function(genes_entrez, bg_entrez, ct_name, outdir,
                             org_db, kegg_org, pval_cut, qval_cut, show_cat) {

  # ---- GO 富集 (BP / CC / MF) ----
  go_bp <- enrichGO(
    gene          = genes_entrez,
    universe      = bg_entrez,
    OrgDb         = org_db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cut,
    qvalueCutoff  = qval_cut,
    readable      = TRUE
  )

  go_cc <- enrichGO(
    gene          = genes_entrez,
    universe      = bg_entrez,
    OrgDb         = org_db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cut,
    qvalueCutoff  = qval_cut,
    readable      = TRUE
  )

  go_mf <- enrichGO(
    gene          = genes_entrez,
    universe      = bg_entrez,
    OrgDb         = org_db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cut,
    qvalueCutoff  = qval_cut,
    readable      = TRUE
  )

  # ---- 汇总 GO ----
  go_all <- NULL
  for (ont_name in c("BP", "CC", "MF")) {
    go_res <- switch(ont_name, BP = go_bp, CC = go_cc, MF = go_mf)
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      df <- as.data.frame(go_res)
      df$ONTOLOGY <- ont_name
      go_all <- rbind(go_all, df)
    }
  }

  if (!is.null(go_all) && nrow(go_all) > 0) {
    # 输出表格
    go_csv <- file.path(outdir, paste0(ct_name, "_GO_enrichment.csv"))
    write_csv(go_all, go_csv)
    message(sprintf("  ✓ GO: %d 条显著富集 → %s", nrow(go_all), go_csv))

    # 气泡图：按 ONTOLOGY 分面展示 top N
    go_top <- go_all %>%
      group_by(ONTOLOGY) %>%
      slice_min(order_by = p.adjust, n = show_cat, with_ties = FALSE) %>%
      ungroup()

    p_go <- ggplot(go_top, aes(x = GeneRatio, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      labs(title = paste0("GO Enrichment — ", ct_name),
           x = "GeneRatio", y = "") +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(face = "bold"))

    go_pdf <- file.path(outdir, paste0(ct_name, "_GO_dotplot.pdf"))
    n_panels <- length(unique(go_top$ONTOLOGY))
    ggsave(go_pdf, p_go, width = 10, height = max(6, n_panels * 3))
    message(sprintf("  ✓ GO dotplot → %s", go_pdf))
  } else {
    message(sprintf("  ⚠ GO: 无显著富集条目"))
    # 写入空表格，确保 Snakemake output 确定性
    go_csv <- file.path(outdir, paste0(ct_name, "_GO_enrichment.csv"))
    write_csv(data.frame(ONTOLOGY=character(), Description=character(),
                         GeneRatio=character(), p.adjust=numeric(),
                         stringsAsFactors=FALSE), go_csv)
    go_pdf <- file.path(outdir, paste0(ct_name, "_GO_dotplot.pdf"))
    pdf(go_pdf, width = 8, height = 4)
    plot.new(); text(0.5, 0.5, "No significant GO enrichment", cex = 1.2)
    dev.off()
    message(sprintf("  ✓ GO placeholder → %s", go_pdf))
  }

  # ---- KEGG 富集 ----
  kegg_res <- enrichKEGG(
    gene          = genes_entrez,
    organism      = kegg_org,
    universe      = bg_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cut,
    qvalueCutoff  = qval_cut
  )

  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_df <- as.data.frame(kegg_res)
    kegg_csv <- file.path(outdir, paste0(ct_name, "_KEGG_enrichment.csv"))
    write_csv(kegg_df, kegg_csv)
    message(sprintf("  ✓ KEGG: %d 条显著富集 → %s", nrow(kegg_df), kegg_csv))

    # 气泡图
    kegg_top <- kegg_df %>%
      slice_min(order_by = p.adjust, n = show_cat, with_ties = FALSE)

    p_kegg <- ggplot(kegg_top, aes(x = GeneRatio, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
      labs(title = paste0("KEGG Pathway Enrichment — ", ct_name),
           x = "GeneRatio", y = "") +
      theme_bw(base_size = 10)

    kegg_pdf <- file.path(outdir, paste0(ct_name, "_KEGG_dotplot.pdf"))
    ggsave(kegg_pdf, p_kegg, width = 10, height = max(5, nrow(kegg_top) * 0.35))
    message(sprintf("  ✓ KEGG dotplot → %s", kegg_pdf))
  } else {
    message(sprintf("  ⚠ KEGG: 无显著富集通路"))
    # 写入空表格，确保 Snakemake output 确定性
    kegg_csv <- file.path(outdir, paste0(ct_name, "_KEGG_enrichment.csv"))
    write_csv(data.frame(Description=character(), GeneRatio=character(),
                         p.adjust=numeric(), stringsAsFactors=FALSE),
              kegg_csv)
    kegg_pdf <- file.path(outdir, paste0(ct_name, "_KEGG_dotplot.pdf"))
    pdf(kegg_pdf, width = 8, height = 4)
    plot.new(); text(0.5, 0.5, "No significant KEGG enrichment", cex = 1.2)
    dev.off()
    message(sprintf("  ✓ KEGG placeholder → %s", kegg_pdf))
  }

  invisible(NULL)
}

# ============================================================
# 主流程
# ============================================================

# 收集所有基因作为 background（用于 universe 参数）
all_genes_entrez <- c()

for (i in seq_along(sig_files)) {
  ct_name  <- contrast_names[i]
  sig_file <- sig_files[i]

  message(sprintf("\n── %s ──────────────────────────────", ct_name))

  if (!file.exists(sig_file)) {
    warning(sprintf("文件不存在，跳过: %s", sig_file))
    next
  }

  sig <- read_csv(sig_file, show_col_types = FALSE)

  if (!(gene_id_col %in% colnames(sig))) {
    stop(sprintf("列 '%s' 不存在于 %s。可用列: %s",
                 gene_id_col, sig_file, paste(colnames(sig), collapse = ", ")))
  }

  gene_ids <- unique(sig[[gene_id_col]])
  message(sprintf("  差异基因: %d 个", length(gene_ids)))

  # 基因 ID → ENTREZID 转换
  entrez_res <- bitr(
    gene_ids,
    fromType = from_type,
    toType   = "ENTREZID",
    OrgDb    = org_db
  )
  genes_entrez <- unique(entrez_res$ENTREZID)
  all_genes_entrez <- unique(c(all_genes_entrez, genes_entrez))
  message(sprintf("  %s→ENTREZID: %d/%d 成功映射", from_type, length(genes_entrez), length(gene_ids)))

  if (length(genes_entrez) == 0) {
    warning(sprintf("  无基因可映射到 ENTREZID，跳过 %s", ct_name))
    next
  }

  enrich_and_plot(
    genes_entrez = genes_entrez,
    bg_entrez    = NULL,   # 先用 NULL，后续可用 all_genes_entrez
    ct_name      = ct_name,
    outdir       = outdir,
    org_db       = org_db,
    kegg_org     = kegg_org,
    pval_cut     = pval_cut,
    qval_cut     = qval_cut,
    show_cat     = show_cat
  )
}

message(sprintf("\n[clusterProfiler] 完成。输出目录: %s", outdir))
