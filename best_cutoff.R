#


# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")



library(survival)
library(survminer)
library(pROC)


# 对于1002批次而言




# 确定连续变量的最佳cutoff值
surv_cutpoint(
  data = d,
  time = time,
  event = status,
  variables = ''
)
surv_categorize()



