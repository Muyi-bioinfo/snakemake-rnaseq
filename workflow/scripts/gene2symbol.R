#!/usr/bin/env Rscript
# ============================================================
# gene2symbol.R — 基因 ID → Symbol 注释工具
#
# 依赖: bioconductor-rtracklayer (GTF 标准解析), r-dplyr, r-readr
#
# 支持两种使用模式:
#
# 1. 独立运行（为任意 CSV 追加 gene symbol 列）:
#    Rscript gene2symbol.R --input DEG.csv --gtf genes.gtf --output DEG_anno.csv
#
# 2. 被 deseq2.R source() 引用:
#    source("workflow/scripts/gene2symbol.R")
#    anno <- load_gtf_annotation("genes.gtf")
#    df   <- add_gene_symbol(df, anno, id_col = "gene_id")
#
# 输入要求:
#   - 输入 CSV 需含 id_col 指定的列（默认 "gene_id"）
#   - GTF 的 gene 行中需同时含 gene_id 和 gene_name 属性
# ============================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(readr)
  library(dplyr)
})


# ============================================================
# 核心函数（source() 后可直接调用）
# ============================================================

#' 从 GTF 文件中提取 gene_id → gene_name 映射
#'
#' 使用 rtracklayer::import() 标准化解析 GTF，自动处理属性格式差异。
#' 仅保留 type=="gene" 的行以提高提取速度。
#'
#' @param gtf_path GTF 文件路径
#' @return data.frame, 列: gene_id, gene_name（均去重）
load_gtf_annotation <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF 文件不存在: ", gtf_path)
  }

  # 只读 gene 行（跳过 exon/CDS/UTR 等，大幅减少内存）
  gtf <- rtracklayer::import(gtf_path, format = "gtf",
                             features = "gene")

  # 提取 gene_id 和 gene_name（mcols 自动解析属性列）
  meta <- as.data.frame(mcols(gtf))

  has_id   <- "gene_id"   %in% colnames(meta)
  has_name <- "gene_name" %in% colnames(meta)

  if (!has_id || !has_name) {
    present <- colnames(meta)
    stop("GTF 中未同时找到 gene_id + gene_name。现有属性: ",
         paste(present, collapse = ", "))
  }

  anno <- data.frame(
    gene_id   = as.character(meta$gene_id),
    gene_name = as.character(meta$gene_name),
    stringsAsFactors = FALSE
  )

  # 去重：每个 gene_id 保留第一个匹配的 gene_name
  anno <- anno[!duplicated(anno$gene_id), ]

  message(sprintf("[gene2symbol] rtracklayer::import → %d 个 gene_id→gene_name 映射",
                  nrow(anno)))

  return(anno)
}


#' 为数据框追加 gene symbol 列
#'
#' 以左连接方式将 gene_id 映射为 gene_name。原 df 中无匹配的行保留原样
#' （gene_name 填入 NA）。
#'
#' @param df      输入数据框，须含 id_col 指定的列
#' @param anno    由 load_gtf_annotation() 返回的注释表
#' @param id_col  df 中基因 ID 的列名 (默认 "gene_id")
#' @return 新 data.frame，在最右追加 gene_name 列
add_gene_symbol <- function(df, anno, id_col = "gene_id") {
  if (!id_col %in% colnames(df)) {
    warning(sprintf("[gene2symbol] 列 '%s' 不存在，跳过注释", id_col))
    return(df)
  }

  n_before <- nrow(df)
  df <- df %>%
    left_join(anno, by = setNames("gene_id", id_col))

  n_matched <- sum(!is.na(df$gene_name))
  message(sprintf("[gene2symbol] 添加 gene_name: %d/%d (%.1f%%) 命中",
                  n_matched, n_before, 100 * n_matched / n_before))

  return(df)
}


# ============================================================
# 独立运行模式
# ============================================================

if (sys.nframe() == 0L) {          # 仅在直接运行 Rscript 时执行

  args <- commandArgs(trailingOnly = TRUE)

  usage <- function() {
    cat("\nUsage:\n")
    cat("  独立运行: Rscript gene2symbol.R --input <csv> --gtf <gtf> [--output <csv>] [--id_col <col>]\n")
    cat("  库引用:   source('gene2symbol.R'); anno <- load_gtf_annotation(...)\n\n")
    quit(status = 1)
  }

  if (length(args) == 0) usage()

  source("workflow/scripts/common.R")

  input_file  <- parse_arg(args, "--input")
  gtf_file    <- parse_arg(args, "--gtf")
  output_file <- parse_arg(args, "--output")
  id_col      <- parse_arg(args, "--id_col")

  if (is.null(input_file) || is.null(gtf_file)) {
    cat("[ERROR] --input 和 --gtf 为必需参数\n")
    usage()
  }
  if (is.null(id_col)) id_col <- "gene_id"
  if (!file.exists(input_file)) stop("输入文件不存在: ", input_file)

  message("========================================")
  message("gene2symbol — 基因 ID → Symbol 注释")
  message("========================================")
  message("输入:  ", input_file)
  message("GTF:   ", gtf_file)
  message("ID 列: ", id_col)
  message("")

  anno <- load_gtf_annotation(gtf_file)

  df <- read_csv(input_file, show_col_types = FALSE, progress = FALSE)
  message(sprintf("[gene2symbol] 读入数据: %d 行 × %d 列", nrow(df), ncol(df)))

  df <- add_gene_symbol(df, anno, id_col = id_col)

  if (is.null(output_file)) {
    output_file <- sub("\\.csv$", "_anno.csv", input_file)
  }
  write_csv(df, output_file)
  message(sprintf("[gene2symbol] 输出:  %s (%d 行 × %d 列)", output_file, nrow(df), ncol(df)))
  message("========================================")
}
