# compare_rf.R
library(TreeDist)
library(ape)

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare_rf.R tree1.nwk tree2.nwk")
}

tree1 <- read.tree(args[1])
tree2 <- read.tree(args[2])
print(tree1)
print(tree2)
print(find("RobinsonFoulds"))

# 计算 RF 距离
rf <- TreeDist::RobinsonFoulds(tree1, tree2)
print(rf)
print(typeof(rf))

if (is.list(rf)) {
  # 新版 TreeDist 返回 list，正常解析
  cat("RF distance (raw):", rf$RawDistance, "\n")
  cat("RF distance (normalized):", rf$NormalizedDistance, "\n")
} else if (is.numeric(rf)) {
  # 旧版或别的函数只返回单个值
  cat("RF distance:", rf, "\n")
} else {
  stop("Unknown RF return type: ", typeof(rf))
}