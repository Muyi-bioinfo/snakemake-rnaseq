# ============================================================
# RNA-seq 流程 — R 脚本公共函数
# 被 deseq2.R / clusterprofiler.R / gene2symbol.R source()
# ============================================================

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx == length(args)) {
    stop(sprintf("参数 %s 缺少值（不能是命令行最后一个参数）", flag))
  }
  args[idx + 1]
}
