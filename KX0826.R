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
library(easystats)
library(geepack)
library(lme4)
library(lmerTest)
library(ggeffects)

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
dose <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002_Medical Review Listing.xlsx',
  sheet = '试验药物_对照药品给药'
)
pk_kx826 <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002-群体暴露水平血液样本采集-20210819.xlsx',
  sheet = 'KX-826 PK Data'
)

pk_kx982 <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/KX0826-CN-1002-群体暴露水平血液样本采集-20210819.xlsx',
  sheet = 'KX-982 PK Data'
)
cls.dose <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/Subjects_Rave_RTSM_08-18-2021_10 23 38.xlsx', skip = 2) %>%
  mutate(`受试者ID`=str_c('0', `受试者ID`))


print(length(unique(drug_his$受试者号))) # 共166例患者
glue::glue('total drugs {length(table(drug_his$药物名称))}') # 竟然204种药

idrugs <- c('米诺地尔酊', '非那雄胺片', '启悦、蓝乐非那雄胺片、仙琚非那雄胺片')
feina_id <- drug_his %>% filter(`药物名称` %in% c('非那雄胺片', '非那雄胺')) %>% pull(`受试者号`)
minuo_id <- drug_his %>% filter(`药物名称` %in% c('米诺地尔酊', '米诺地尔')) %>% pull(`受试者号`)
intersect(minuo_id, feina_id)


foo <- clin_diag_studys %>% select(`受试者号`, `访视名称`, `目标区域内非毳毛数（TAHC）`) %>%
  pivot_wider(names_from = `访视名称`, values_from = `目标区域内非毳毛数（TAHC）`)


hga %>%
  mutate(`头发生长情况（HGA）评估-研究者评估`=str_extract(`头发生长情况（HGA）评估-研究者评估`, '-?\\d+')) %>%
  select(`受试者号`,`访视名称`, `头发生长情况（HGA）评估-研究者评估`)

jisu2 <- jisu %>% select(`受试者号`,`血清睾酮`, `访视名称`) %>%
  pivot_wider(names_from = `访视名称`, values_from = `血清睾酮`, names_prefix = 'jisu_'
  )
pk_kx826.2 <- pk_kx826 %>% select(`受试者编号`, `时间点`, `检测结果`) %>%
  pivot_wider(names_from = `时间点`, values_from = `检测结果`, names_prefix = 'kx826_')
pk_kx982.2 <- pk_kx982 %>% select(`受试者编号`, `时间点`, `检测结果`) %>%
  pivot_wider(names_from = `时间点`, values_from = `检测结果`, names_prefix = 'kx928_')
cls.dose2 <- cls.dose %>% select(`受试者ID`, `研究项目分组`)

bar <- patient_summary %>% select(`受试者号`, `年龄`, `性别`, `BMI`) %>%
  left_join(foo, by = c('受试者号')) %>%
  left_join((clin_diag_study0 %>% select(`受试者号`, `雄激素性秃发严重程度检查结果`)),
            by = c('受试者号')) %>%
  left_join((jisu2), by = c('受试者号')) %>%
  left_join(pk_kx826.2, by = c('受试者号'='受试者编号')) %>%
  left_join(pk_kx982.2, by = c('受试者号'='受试者编号')) %>%
  left_join(cls.dose2, by = c('受试者号'='受试者ID'))


# write_excel_csv(bar, '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/meibo_PK.csv')

# 以毛发检验为outcome的分析 --------------------------------------------------------

# 针对foo的重复测量数据分别进行gee和lme, foo为具有结局指标的120例样本

foo2 <- foo %>% left_join((patient_summary %>% select(`受试者号`, `年龄`, `BMI`)),
                          by = c('受试者号')) %>% select(-`计划外访视`) %>%
  left_join(cls.dose2, by = c('受试者号'='受试者ID')) %>%
  left_join((clin_diag_study0 %>% select(`受试者号`, `雄激素性秃发严重程度检查结果`)),
            by = c('受试者号')) %>%
  rename('dislevel'='雄激素性秃发严重程度检查结果')

df <- foo2 %>% pivot_longer(
  cols = -c(`受试者号`, `年龄`, BMI, `研究项目分组`, W1D1, dislevel),
  names_to = 'time',
  values_to = 'score'
) %>%
  rename('PatientID'='受试者号',
         'age' = '年龄',
         'group' = '研究项目分组'
         ) %>%
  mutate(across(c(score, W1D1, BMI, age), as.numeric))
df$time <- factor(df$time, levels = c('W6D1','W12D1', 'W18D1', 'W24D1'))
df$group <- factor(df$group, levels = c('安慰剂QD', '安慰剂BID','2.5mgBID','5mgQD','5mgBID'))


# gee only 主效应
geefit_main <- geeglm(
  score ~ time + group,
  id=PatientID,
  corstr='exchangeable',
  family="gaussian",
  data=df,
  std.err = 'san.se'
)
summary(geefit_main)
anova(geefit_main)
broom::tidy(geefit_main)

# gee 考虑所有的协变量
geefit2 <- geeglm(score ~ time + group + age + BMI + W1D1 + dislevel,
                  id=PatientID,
                  corstr='exchangeable',
                  family="gaussian",
                  data=df,
                  std.err = 'san.se')
summary(geefit2)
QIC(geefit2)

anova(geefit_main, geefit2)

## 多水平模型/线性混合模型
lme_model <-
  lmerTest::lmer(score ~ time + group + age + BMI + W1D1 + dislevel + (1|`PatientID`),
             data = df
  )

summary(lme_model)
lmerTest::ranova(lme_model)
anova(lme_model, type = 'III', ddf="Satterthwaite")
aa <- ggpredict(lme_model, 'time')

plot(ggemmeans(lme_model, terms = c("time", "group"),
               condition = c(diagnose = "severe")),
     # facets = T
     ) +
  ggplot2::ggtitle("GLMER Effect plot")

##### plot model by estimate
# sjPlot::plot_model(lme_model)
plot(parameters(lme_model))


# 不同的分组对比 -----------------------------------------------------------------

# 只选择BID组的数据
df2 <- df %>% filter(group %in% c('安慰剂BID','2.5mgBID','5mgBID'))
df2$group <-  factor(df2$group, levels = c('安慰剂BID','2.5mgBID','5mgBID'))

# only 主效应
geefit_main2 <- geeglm(
  score ~ time + group,
  id=PatientID,
  corstr='exchangeable',
  family="gaussian",
  data=df2,
  std.err = 'san.se'
)
summary(geefit_main2)
anova(geefit_main2)

# 考虑所有的协变量
geefit3 <- geeglm(score ~ time + group + age + BMI + W1D1 + dislevel,
                  id=PatientID,
                  corstr='exchangeable',
                  family="gaussian",
                  data=df2,
                  std.err = 'san.se')
summary(geefit3)
anova(geefit3)
QIC(geefit3)

anova(geefit_main2, geefit3)


# lme model
lme_model2 <-
  lmerTest::lmer(score ~ time + group + age + BMI + W1D1 + dislevel + (1|`PatientID`),
             data = df2
  )

summary(lme_model2)
anova(lme_model2)

plot(ggemmeans(lme_model2, terms = c("time", "group"),
               condition = c(diagnose = "severe"))) +
  ggplot2::ggtitle("GLMER Effect plot")


# normal statistic for group BID ------------------------------------------





