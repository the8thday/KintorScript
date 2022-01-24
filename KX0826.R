## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2022-01-24
##
## Copyright (c) cliu, 2022
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 此为KX0826-CN-1002 的患者入组数据。
##
##
## ---------------------------



library(tidyverse)


# tidy data ---------------------------------------------------------------


drug_his <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '既往_合并用药'
)
clin_diag_study0 <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '雄激素性秃发临床诊断及严重程度'
)
clin_diag_studys <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '毛发检查'
  )
hga <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '毛发生长情况评估_HGA'
  ) # 此表没问题吧？
hga_third <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = 'HGA_第三方专业医生评估'
)
group_blood <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '群体暴露水平血液样本采集'
)
patient_summary <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '人口统计学信息'
  )
jisu <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '激素检查'
)

print(length(unique(drug_his$受试者号))) # 共166例患者
glue::glue('total drugs {length(table(drug_his$药物名称))}') # 竟然204种药

idrugs <- c('米诺地尔酊', '非那雄胺片', '启悦、蓝乐非那雄胺片、仙琚非那雄胺片')
drug_his %>% filter(`药物名称` %in% c('非那雄胺片', '非那雄胺'))



# 以毛发检验为outcome的分析 --------------------------------------------------------

foo <- clin_diag_studys %>% select(`受试者号`, `访视名称`, `目标区域内非毳毛数（TAHC）`) %>%
  pivot_wider(names_from = `访视名称`, values_from = `目标区域内非毳毛数（TAHC）`)


hga %>%
  mutate(`头发生长情况（HGA）评估-研究者评估`=str_extract(`头发生长情况（HGA）评估-研究者评估`, '-?\\d+')) %>%
  select(`受试者号`,`访视名称`, `头发生长情况（HGA）评估-研究者评估`)


bar <- patient_summary %>% select(`受试者号`, `年龄`, `性别`, `BMI`) %>%
  left_join(foo, by = c('受试者号')) %>%
  left_join((clin_diag_study0 %>% select(`受试者号`, `雄激素性秃发严重程度检查结果`)),
            by = c('受试者号')) %>%
  left_join((jisu %>% select(`受试者号`,`血清睾酮`)), by = c('受试者号'))

# write_delim(bar, '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/meibo_PK.txt', delim = '\t')







