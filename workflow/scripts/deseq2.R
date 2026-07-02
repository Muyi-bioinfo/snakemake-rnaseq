#!/usr/bin/env Rscript
# ============================================================
# DESeq2 差异表达分析脚本
# 输入：featureCounts 计数矩阵 + 样本分组信息
# 输出：差异基因列表 + 诊断图表
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

set.seed(42)  # 确保 DESeq2 结果可复现

# ---- 解析命令行参数 ----
source("workflow/scripts/common.R")

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("Usage: Rscript deseq2.R --counts <file> --groups <json> --contrasts <json> \\
       --outdir <dir> [--gtf <gtf>] --padj <float> --log2fc <float> --top_n <int>\n")
  quit(status = 1)
}

if (length(args) == 0) usage()

counts_file  <- parse_arg(args, "--counts")
groups_json  <- parse_arg(args, "--groups")
contr_json   <- parse_arg(args, "--contrasts")
outdir       <- parse_arg(args, "--outdir")
gtf_file     <- parse_arg(args, "--gtf")
padj_cut     <- as.numeric(parse_arg(args, "--padj"))
log2fc_cut   <- as.numeric(parse_arg(args, "--log2fc"))
top_n        <- as.integer(parse_arg(args, "--top_n"))

# 默认值
if (is.null(padj_cut))   padj_cut   <- 0.05
if (is.null(log2fc_cut)) log2fc_cut <- 1.0
if (is.null(top_n))      top_n      <- 50
if (is.null(outdir))     outdir     <- "results/06_deseq2"

# gene symbol 注释（可选）
has_gtf <- !is.null(gtf_file) && file.exists(gtf_file)
if (has_gtf) {
  source("workflow/scripts/gene2symbol.R")
  anno <- load_gtf_annotation(gtf_file)
  message("")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[DESeq2] ========================================")
message("[DESeq2] DESeq2 差异表达分析")
message("[DESeq2] ========================================")
message("[DESeq2] Counts file: ", counts_file)
message("[DESeq2] Output dir:  ", outdir)
message("[DESeq2] padj cutoff: ", padj_cut)
message("[DESeq2] log2FC cutoff: ", log2fc_cut)
message("[DESeq2]")

# ---- 1. 读入计数矩阵 ----
message("[DESeq2] [1/6] 读入 featureCounts 计数矩阵...")
counts_raw <- read.table(counts_file, header = TRUE, row.names = 1,
                         comment.char = "#", check.names = FALSE)

# featureCounts 输出格式：前几列是基因注释信息(Geneid, Chr, Start, End, Strand, Length)
# 后面每列是一个样本的计数
count_cols <- grep("\\.sorted\\.bam$", colnames(counts_raw))
if (length(count_cols) > 0) {
  counts <- counts_raw[, count_cols, drop = FALSE]
  gene_info <- counts_raw[, -count_cols, drop = FALSE]
} else {
  # fallback: 假设第7列以后是counts
  counts <- counts_raw[, 7:ncol(counts_raw), drop = FALSE]
  gene_info <- counts_raw[, 1:6, drop = FALSE]
}

# 简化样本名
colnames(counts) <- gsub("\\.sorted\\.bam$", "", colnames(counts))
message(sprintf("[DESeq2]   基因数: %d, 样本数: %d", nrow(counts), ncol(counts)))

# ---- 2. 构建样本元数据 ----
message("[DESeq2] [2/6] 构建样本元数据...")
groups <- fromJSON(groups_json)   # e.g. {"WT":["WT_1","WT_2","WT_3"], "KO":["KO_1","KO_2","KO_3"]}

# 展开为 sample → group 映射（向量化，避免 rbind 逐次拷贝）
sample_group <- data.frame(
  sample = unlist(groups, use.names = FALSE),
  group  = rep(names(groups), lengths(groups)),
  stringsAsFactors = FALSE
)

# 只保留在计数矩阵中存在的样本
sample_group <- sample_group[sample_group$sample %in% colnames(counts), ]
rownames(sample_group) <- sample_group$sample

# 按计数矩阵列顺序重排
sample_group <- sample_group[colnames(counts), , drop = FALSE]

message(sprintf("[DESeq2]   分组: %s", paste(unique(sample_group$group), collapse = ", ")))
print(sample_group)

# ---- 3. 构建 DESeqDataSet ----
message("[DESeq2] [3/6] 构建 DESeqDataSet...")
sample_group$group <- factor(sample_group$group)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = sample_group,
  design    = ~ group
)

# 预过滤：去除低表达基因
keep <- rowSums(counts(dds) >= 10) >= min(table(sample_group$group))
dds <- dds[keep, ]
gene_info <- gene_info[keep, , drop = FALSE]   # 同步过滤注释表，保持行对齐
message(sprintf("[DESeq2]   过滤后基因数: %d (过滤前: %d)", nrow(dds), nrow(counts)))

# ---- 4. 运行 DESeq2 ----
message("[DESeq2] [4/6] 运行 DESeq2 差异分析...")
dds <- DESeq(dds)

# VST 变换（用于可视化的标准化表达量）
vsd <- vst(dds, blind = TRUE)
message("[DESeq2]   DESeq2 分析完成")

# ---- 5. 提取差异基因 ----
message("[DESeq2] [5/6] 提取差异基因...")
contrasts <- fromJSON(contr_json)
# e.g. [{"name":"KO_vs_WT","case":"KO","control":"WT"}]

all_degs <- list()

for (i in seq_len(nrow(contrasts))) {
  ct_name    <- contrasts$name[i]
  ct_case    <- contrasts$case[i]
  ct_control <- contrasts$control[i]

  message(sprintf("[DESeq2]   处理对比: %s (%s vs %s)", ct_name, ct_case, ct_control))

  res <- results(dds, contrast = c("group", ct_case, ct_control),
                 alpha = padj_cut, independentFiltering = TRUE)
  res <- res[order(res$pvalue), ]

  # 添加基因注释信息（按 gene_id merge，避免位置错位）
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  gene_info$gene_id <- rownames(gene_info)
  res_df <- merge(res_df, gene_info, by = "gene_id", all.x = TRUE, sort = FALSE)

  # 筛选显著性基因
  sig <- na.omit(res_df[res_df$padj < padj_cut &
                         abs(res_df$log2FoldChange) >= log2fc_cut, ])
  sig <- sig[order(sig$padj), ]

  # ---- gene symbol 注释（可选）----
  if (has_gtf) {
    res_df <- add_gene_symbol(res_df, anno)
    sig    <- add_gene_symbol(sig,    anno)
  }

  # 输出完整结果表
  out_file <- file.path(outdir, paste0(ct_name, "_all_results.csv"))
  write_csv(res_df, out_file)

  sig_file <- file.path(outdir, paste0(ct_name, "_significant.csv"))
  write_csv(sig, sig_file)

  # 统计
  n_up   <- sum(sig$log2FoldChange > 0)
  n_down <- sum(sig$log2FoldChange < 0)
  message(sprintf("[DESeq2]     DEGs: %d (↑%d ↓%d)", nrow(sig), n_up, n_down))

  all_degs[[ct_name]] <- list(
    res = res,
    sig = sig,
    case = ct_case,
    control = ct_control
  )
}

# ---- 6. 诊断图表 ----
message("[DESeq2] [6/6] 生成诊断图表...")

# --- 6a. PCA 图 ---
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title = "PCA — 样本间表达量距离",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(outdir, "PCA_plot.pdf"), p_pca, width = 7, height = 6)
ggsave(file.path(outdir, "PCA_plot.png"), p_pca, width = 7, height = 6, dpi = 150)

# --- 6b. 样本距离热图 ---
sample_dist <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)
rownames(sample_dist_mat) <- colnames(sample_dist_mat) <- colnames(vsd)

pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         main = "样本间表达量距离 (VST)",
         display_numbers = TRUE, number_format = "%.0f")
dev.off()

# --- 6c. 每个对比的 MA plot + Volcano plot ---
for (ct_name in names(all_degs)) {
  info <- all_degs[[ct_name]]
  res_df <- as.data.frame(info$res)
  res_df$gene_id <- rownames(res_df)
  res_df$sig <- ifelse(
    !is.na(res_df$padj) & res_df$padj < padj_cut &
    abs(res_df$log2FoldChange) >= log2fc_cut,
    ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "NS"
  )
  res_df$label <- ifelse(res_df$sig != "NS", res_df$gene_id, NA)

  # MA plot
  p_ma <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = sig)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
    geom_hline(yintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", alpha = 0.4) +
    labs(title = paste("MA plot —", ct_name),
         x = "Mean of normalized counts", y = "log2 Fold Change") +
    theme_bw(base_size = 13)

  ggsave(file.path(outdir, paste0(ct_name, "_MA_plot.pdf")), p_ma, width = 7, height = 5)

  # Volcano plot
  res_df$neglog10_padj <- -log10(res_df$padj)
  top_labels <- res_df[order(res_df$padj), ]
  top_labels <- head(top_labels[top_labels$sig != "NS", ], 15)
  # 优先使用 gene_name 标签
  if (!is.null(top_labels$gene_name)) {
    top_labels$label <- ifelse(is.na(top_labels$gene_name),
                               top_labels$gene_id, top_labels$gene_name)
  } else {
    top_labels$label <- top_labels$gene_id
  }

  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj, color = sig)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
    geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", alpha = 0.4) +
    geom_text_repel(data = top_labels, aes(label = label), size = 2.5,
                    max.overlaps = 20, show.legend = FALSE) +
    labs(title = paste("Volcano plot —", ct_name),
         x = "log2 Fold Change", y = "-log10 adjusted p-value") +
    theme_bw(base_size = 13)

  ggsave(file.path(outdir, paste0(ct_name, "_volcano_plot.pdf")), p_volcano, width = 7, height = 6)
}

# --- 6d. 显著差异基因热图 (top N) ---
if (length(all_degs) > 0) {
  # 合并所有对比中的显著基因
  all_sig_genes <- unique(unlist(lapply(all_degs, function(x) {
    rownames(x$sig)[1:min(top_n, nrow(x$sig))]
  })))
  all_sig_genes <- all_sig_genes[!is.na(all_sig_genes)]

  if (length(all_sig_genes) > 2 && length(all_sig_genes) <= 200) {
    mat <- assay(vsd)[all_sig_genes, , drop = FALSE]
    mat <- mat - rowMeans(mat)  # centered

    pdf(file.path(outdir, "DEG_heatmap.pdf"), width = max(8, ncol(vsd) * 1.2),
        height = max(6, min(length(all_sig_genes) * 0.15, 15)))
    pheatmap(mat,
             annotation_col = data.frame(
               group = colData(dds)$group,
               row.names = colnames(dds)
             ),
             show_rownames = length(all_sig_genes) <= 60,
             main = paste("Top DEGs Heatmap (n =", length(all_sig_genes), ")"))
    dev.off()
    message(sprintf("[DESeq2]   热图基因数: %d", length(all_sig_genes)))
  }
}

# ---- 输出汇总 ----
message("[DESeq2]")
message("[DESeq2] ========================================")
message("[DESeq2] DESeq2 分析完成！")
message("[DESeq2] 输出目录: ", outdir)
message("[DESeq2] 输出文件:")
for (f in list.files(outdir, pattern = "\\.(csv|pdf|png)$")) {
  message(sprintf("[DESeq2]   - %s", f))
}
message("[DESeq2] ========================================")
