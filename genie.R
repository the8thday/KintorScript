library(GENIE3)


# 用法 ----------------------------------------------------------------------

GENIE3(
  exprMatrix, # 表达矩阵
  regulators = NULL, # 指定潜在的调控因子，比如转录因子等
  targets = NULL, # 潜在的被调控的靶标基因
  treeMethod = "RF", # 选择方法，默认的是“RF"（随机森林），还可以选择“ET”（Extra-Trees）
  K = "sqrt", 
  nTrees = 1000, # 树的量，默认是1000
  nCores = 1, # 用于并行计算的核数，表达矩阵较大时选择并行，运算速度更快。
  returnMatrix = TRUE, # 结果返回形式是矩阵还是list，选择"TRUE"就返回矩阵，否则就返回list
  verbose = FALSE # 是否展示计算进度，默认是FALSE，即不展示计算进度
)





