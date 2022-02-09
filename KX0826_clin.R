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
## Notes: 此为KX0826-CN-1002 的患者入组数据。主要采用的分析方法为重复测量资料的数据分析
##
##
## ---------------------------

library(tidyverse)
library(easystats)
library(geepack)
library(lme4)
# library(brms)
library(lmerTest)
# library(merTools)
# library(mixedup) # 一个提取mixed model的小工具
library(ggeffects)
library(emmeans)

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

idrugs <- c('米诺地尔酊', '非那雄胺片', '启悦','蓝乐', '仙琚')
feina_id <- drug_his %>% filter(str_detect(`药物名称`, '非那雄胺')) %>%
  pull(`受试者号`)
minuo_id <- drug_his %>% filter(str_detect(`药物名称`, '米诺地尔')) %>% pull(`受试者号`)
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
  rename('dislevel'='雄激素性秃发严重程度检查结果') %>%
  mutate(feina = ifelse(`受试者号` %in% feina_id, 1, 0),
         minuo = ifelse(`受试者号` %in% minuo_id, 2, 0),
         drughis = feina + minuo
         ) %>% select(-c(feina, minuo)) %>%
  mutate(across(starts_with('W'), as.numeric),
         diffW6 = W6D1 - W1D1,
         diffW12 = W12D1 - W1D1,
         diffW18 = W18D1 - W1D1,
         diffW24 = W24D1 - W1D1
         ) %>%
  filter(!`受试者号` %in% c('01004', '03032', '03034', '06007'))

df_p <- foo2 %>% pivot_longer(
  cols = starts_with('W'),
  names_to = 'time',
  values_to = 'score'
) %>%
  rename('PatientID'='受试者号',
         'age' = '年龄',
         'group' = '研究项目分组'
         ) %>%
  mutate(across(c(score, BMI, age), as.numeric))
df_p$time <- factor(df_p$time, levels = c('W1D1','W6D1','W12D1', 'W18D1', 'W24D1'))
df_p$group <- factor(df_p$group, levels = c('安慰剂QD', '安慰剂BID','2.5mgBID','5mgQD','5mgBID'))

df <- foo2 %>%
  pivot_longer(
  cols = -c(`受试者号`, `年龄`, BMI, `研究项目分组`, dislevel, drughis, starts_with('W')),
  names_to = 'time',
  values_to = 'score'
) %>%
  rename('PatientID'='受试者号',
         'age' = '年龄',
         'group' = '研究项目分组'
  ) %>%
  mutate(across(c(score, BMI, age), as.numeric)) %>%
  rstatix::convert_as_factor(drughis)
df$time <- factor(df$time, levels = c('diffW6','diffW12','diffW18', 'diffW24'))

df$group <- factor(df$group, levels = c('安慰剂QD', '安慰剂BID','2.5mgBID','5mgQD','5mgBID'))
# write_excel_csv(df2, '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/gee_lme_BID.csv')


# 针对df的数据探索 ---------------------------------------------------------------

library(showtext)
showtext_auto(enable = TRUE)
# font_add_google('')

df_pp <- foo2 %>% mutate(diffW1=W1D1-W1D1) %>%
  pivot_longer(
    cols = -c(`受试者号`, `年龄`, BMI, `研究项目分组`, dislevel, drughis, starts_with('W')),
    names_to = 'time',
    values_to = 'score'
  ) %>%
  rename('PatientID'='受试者号',
         'age' = '年龄',
         'group' = '研究项目分组'
  ) %>%
  mutate(across(c(score, BMI, age), as.numeric)) %>%
  rstatix::convert_as_factor(drughis)
df_pp$time <- factor(df_pp$time, levels = c('diffW1','diffW6','diffW12','diffW18', 'diffW24'))
df_pp$group <- factor(df_pp$group, levels = c('安慰剂QD', '安慰剂BID','2.5mgBID','5mgQD','5mgBID'))

ggplot(df_pp, mapping = aes(x = time, y = score, color = group, group=group)) +
  stat_summary(geom = 'line',
               fun = 'mean'
  ) +
  stat_summary(geom = 'point',
               fun = 'mean'
  ) +
  # stat_summary(geom = 'errorbar',
  #              fun.min = function(x){mean(x)-sd(x)},
  #              fun.max = function(x){mean(x)+sd(x)},
  #              width = 0.2
  # ) +
  # scale_y_continuous(expand = expansion(mult = c(0,0)),
  #                    limits = c(0, NA)
  # ) +
  ggsci::scale_color_lancet() +
  ggprism::theme_prism() +
  theme(text = element_text(family='STHeiti'))


set.seed(42)
ss <- sample(unique(df$PatientID), 10)

ggplot(df_pp %>% filter(PatientID %in% ss),
       aes(x = time, y = score, color = PatientID)) +
  geom_point() +
  geom_line(aes(group = PatientID)) +
  # facet_grid(Subject ~ ., scales = 'free') +
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0, NA)
  ) +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(21)) +
  ggprism::theme_prism() +
  # theme_bw() +
  theme(text = element_text(family='STHeiti'),
        # axis.text.x = element_text(angle = 90, size = 6),
        legend.position="bottom"
        # axis.text.y = element_text(size = 4),
        # strip.text.y = element_text(angle = 0, size = 6),
        # legend.position = 'none'
  )

ggplot(df_p,
       aes(x = time, y = score, color = PatientID)) +
  geom_point() +
  geom_line(aes(group = PatientID)) +
  facet_grid(group ~ ., scales = 'free') +
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0, NA)
  ) +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(120)) +
  # ggprism::theme_prism() +
  theme_bw() +
  theme(text = element_text(family='STHeiti'),
        # axis.text.x = element_text(angle = 90, size = 6),
        # legend.position="bottom"
        # axis.text.y = element_text(size = 4),
        # strip.text.y = element_text(angle = 0, size = 6),
        legend.position = 'none'
  )

# df 数据的gee lme ----------------------------------------------------------


# 重复测量方差分析
df1 <- df %>% drop_na()
two.way <- rstatix::anova_test(score ~ time * group + Error(`PatientID`/time),
                      data = df1
)
rstatix::get_anova_table(two.way)
pwc <- df1 %>%
  group_by(time) %>%
  pairwise_t_test(
    score ~ group,
    # paired = TRUE,
    p.adjust.method = "fdr"
  )
pwc

# if you don't have a completely balanced User*A*B design, then aov() likely gets in trouble
aov_fit <- aov(score ~ time * group + Error(`PatientID`/time),
           data = df1)
# report::report(fit)
summary(aov_fit)
anova_summary(aov_fit)
parameters::model_parameters(aov_fit)
effectsize::eta_squared(aov_fit)


# gee 考虑所有的协变量的主效应, gee似乎需要na.omit...
geefit_main <- geeglm(
  score ~ time + group + W1D1 + age + BMI + dislevel + drughis,
  id=PatientID,
  corstr='independence',
  family="gaussian",
  data=df,
  std.err = 'san.se'
)
summary(geefit_main)
anova(geefit_main)
broom::tidy(geefit_main)
# Estimated marginal means
plot(ggemmeans(geefit_main, terms = c("time", "group"),
               condition = c(diagnose = "severe")),
     # facets = T
) +
  ggplot2::ggtitle("GEE Effect plot")

# gee 考虑 time score的交互效应
geefit2 <- geeglm(score ~ time * group + W1D1 + age + BMI + dislevel + drughis,
                  id=PatientID,
                  corstr='independence',
                  family="gaussian",
                  data=df,
                  std.err = 'san.se')
summary(geefit2)
anova(geefit2) %>% as.data.frame() %>% rownames_to_column() %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/gee_main.csv')

anova(geefit_main, geefit2) # 结果显示交互效应并没有实际的意义
sjPlot::plot_model(geefit_main)
gtsummary::tbl_regression(geefit_main)

# final gee model
gee_final <- geeglm(score ~ time + group + W1D1 + drughis,
                    id=PatientID,
                    corstr='independence',
                    family="gaussian",
                    data=df,
                    std.err = 'san.se')
summary(gee_final)
QIC(gee_final)
QIC(geefit2)
QIC(geefit_main)
anova(gee_final)
anova(gee_final, geefit_main)
plot(ggemmeans(gee_final, terms = c("time", "group"),
               condition = c(diagnose = "severe")),
     # facets = T
) +
  ggplot2::ggtitle("GEE Effect plot")
gtsummary::tbl_regression(gee_final)


## 多水平模型/线性混合模型 ###
lme_model <-
  lmerTest::lmer(score ~ time + group + W1D1 + drughis +(1|`PatientID`),
             data = df,
             REML = TRUE
  )

summary(lme_model)
confint(lme_model)
head(coef(lme_model)$PatientID)
ranef(coef(lme_model)$PatientID)
performance(lme_model)
lmerTest::ranova(lme_model)
anova(lme_model, type = 'III', ddf="Satterthwaite") # typeIII 模型效应是矫正了其他因素后的各因素主效应结果
car::Anova(lme_model)
# nlme::anova.lme(lme_model, type = "marginal", adjustSigma = F) # 不适用于此函数
aa <- ggpredict(lme_model, 'time')

plot(ggemmeans(lme_model, terms = c("time", "group"),
               condition = c(diagnose = "severe")),
     # facets = T
     ) +
  ggplot2::ggtitle("GLMER Effect plot")

##### plot model by estimate
sjPlot::plot_model(lme_model)
plot(parameters(lme_model))
gtsummary::tbl_regression(lme_model, intercept=TRUE)

# 纳入所有变量
lme_model2 <-
  lmerTest::lmer(score ~ time + group + W1D1 + age + BMI + dislevel + drughis + (1|`PatientID`),
                 data = df
  )

summary(lme_model2)
performance(lme_model2)
lmerTest::ranova(lme_model2)
anova(lme_model2, type = 'III', ddf="Satterthwaite") %>% as.data.frame() %>% rownames_to_column() %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/lme_main.csv')
gtsummary::tbl_regression(lme_model2, intercept=TRUE)


# 不同的分组对比 -----------------------------------------------------------------

# 只选择BID组的数据
df2 <- df %>% filter(group %in% c('安慰剂BID','2.5mgBID','5mgBID'))
df2$group <-  factor(df2$group, levels = c('安慰剂BID','2.5mgBID','5mgBID'))

# only 主效应
geefit_main2 <- geeglm(
  score ~ time + group + W1D1 + drughis,
  id=PatientID,
  corstr='independence',
  family="gaussian",
  data=df2,
  std.err = 'san.se'
)
summary(geefit_main2)
anova(geefit_main2)
gtsummary::tbl_regression(geefit_main2)
anova(geefit_main2, type = 'III', ddf="Satterthwaite") # typeIII 模型效应是矫正了其他因素后的各因素主效应结果
plot(ggemmeans(geefit_main2, terms = c("time", "group"),
               condition = c(diagnose = "severe")),
     # facets = T
) +
  ggplot2::ggtitle("GEE Effect plot")


# lme model
lme_model3 <-
  lmerTest::lmer(score ~ time + group + W1D1 + drughis + (1|`PatientID`),
             data = df2
  )

summary(lme_model3)
anova(lme_model3)
anova(lme_model3, type = 'III', ddf="Satterthwaite")
car::Anova(lme_model3)

plot(ggemmeans(lme_model3, terms = c("time", "group"),
               condition = c(diagnose = "severe"))) +
  ggplot2::ggtitle("GLMER Effect plot")
gtsummary::tbl_regression(lme_model3)

# lme model, 不包含其他变量
lme_model4 <-
  lmerTest::lmer(score ~ time + group + W1D1 + (1|`PatientID`),
                 data = df2
  )

summary(lme_model4)
anova(lme_model4)
anova(lme_model4, type = 'III', ddf="Satterthwaite")
car::Anova(lme_model4)

plot(ggemmeans(lme_model4, terms = c("time", "group"),
               condition = c(diagnose = "severe"))) +
  ggplot2::ggtitle("GLMER Effect plot")
gtsummary::tbl_regression(lme_model4)

# normal statistic for group BID ------------------------------------------

library(rstatix)
# 以6W和24W作为主要研究终点，进行基础统计推断
foo3 <- foo2 %>%
  filter(!`受试者号` %in% c('01004', '03032', '03034', '06007')) %>%
  filter(`研究项目分组` %in% c('安慰剂BID','2.5mgBID','5mgBID')) %>%
  dplyr::select(!starts_with('diff')) %>%
  pivot_longer(
    cols = -c(`受试者号`, `年龄`, BMI, `研究项目分组`, dislevel, drughis),
    names_to = 'time',
    values_to = 'score'
  ) %>%
  fill(score, .direction = 'down') %>%
  pivot_wider(names_from = time, values_from = score) %>%
  mutate(
    across(starts_with('W'), as.numeric),
    diffW6 = W6D1 - W1D1,
    diffW12 = W12D1 - W1D1,
    diffW18 = W18D1 - W1D1,
    diffW24 = W24D1 - W1D1
  )
foo3$研究项目分组 <- factor(foo3$研究项目分组, levels = c('安慰剂BID', '2.5mgBID', '5mgBID'))


# 检验正态性以及方差齐性
foo3 %>% group_by(`研究项目分组`) %>% identify_outliers(diffW24)
foo3 %>% group_by(`研究项目分组`) %>% shapiro_test(diffW24)
foo3 %>% levene_test(diffW24 ~ `研究项目分组`)

# 基础的方差分析
anova_test(data = foo3,
       formula = diffW24 ~ `研究项目分组`
       )

aov.res <- aov(diffW24 ~ `研究项目分组`,
               data = foo3)
report::report(aov.res)
summary(aov.res)
aggregate(foo3$diffW24, by = list(foo3$研究项目分组), FUN = mean)
aggregate(foo3$diffW24, by = list(foo3$研究项目分组), FUN = median)
tuk <- TukeyHSD(aov.res, conf.level = 0.95)
plot(tuk)

(pwc <- foo3 %>%
  pairwise_t_test(
    diffW24 ~ `研究项目分组`,
    paired = FALSE,
    p.adjust.method = "bonferroni"
  ))

(pwc2 <- foo3 %>%
    wilcox_test(
      diffW24 ~ `研究项目分组`,
      paired = FALSE,
      p.adjust.method = "bonferroni"
    ))

# 协方差分析 ANCOVA
ggpubr::ggscatter(
  foo3, x = "W1D1", y = "W24D1",
  color = "研究项目分组", add = "reg.line"
)+
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = `研究项目分组`)
  ) # 并不是很满足线性要求
aov.xie <- aov(diffW24 ~ W1D1 + `研究项目分组`,
               data = foo3)
aov.xie2 <- aov(W24D1 ~ W1D1 + `研究项目分组`,
               data = foo3)
res.aov <- foo3 %>% anova_test(W24D1 ~ W1D1 + `研究项目分组`)
get_anova_table(res.aov) %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/anova.csv')
myfit1 <- lm(diffW24 ~ W1D1 + `研究项目分组`,
             data = foo3)
broom::augment(myfit1)

car::Anova(aov.xie2, type = "III")
car::Anova(myfit1, type = "III")
check_model(aov.xie)
summary(aov.xie)
summary(aov.xie2)
plot(ggemmeans(aov.xie2, terms = c("研究项目分组"),
               condition = c(diagnose = "severe")),
     # facets = T
) +
  ggplot2::ggtitle("ANCOVA Effect plot")
# 两两比较
pwc3 <- emmeans_test(
  W24D1 ~ `研究项目分组`,
  covariate = W1D1,
  p.adjust.method = "fdr",
  data = foo3
)
(pwc3 %>% write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/emmeans_test.csv'))
get_emmeans(pwc3) %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/emmeans_W24.csv')
# 针对W6周的数据
aov.xie3 <- aov(W24D1 ~ W1D1 + `研究项目分组`,
                data = foo3)
summary(aov.xie3)
broom::tidy(aov.xie3) %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/aov_bar.csv')
pwc4 <- emmeans_test(
  W24D1 ~ `研究项目分组`,
  covariate = W1D1,
  p.adjust.method = "fdr",
  data = foo3
)
pwc4 %>% write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/emmeans_bar.csv')
get_emmeans(pwc4) %>%
  write_excel_csv(file = '/Users/congliu/OneDrive/kintor/Daily_Work/KX826/emmeans_foo.csv')


# 只考虑两组
foo4 <- foo3 %>%
  mutate(`研究项目分组` = as.character(`研究项目分组`)) %>%
  filter(`研究项目分组` %in% c('安慰剂BID','5mgBID'))

t_test(data = foo4,
       formula = diffW24 ~ `研究项目分组`,
       paired = F
       )


# 数据填充后的lme分析
df4 <- foo3 %>%
  pivot_longer(
    cols = -c(`受试者号`, `年龄`, BMI, `研究项目分组`, dislevel, drughis, starts_with('W')),
    names_to = 'time',
    values_to = 'score'
  ) %>%
  rename('PatientID'='受试者号',
         'age' = '年龄',
         'group' = '研究项目分组'
  ) %>%
  mutate(across(c(score, BMI, age), as.numeric)) %>%
  rstatix::convert_as_factor(drughis)
df4$time <- factor(df4$time, levels = c('diffW6','diffW12','diffW18', 'diffW24'))
df4$group <- factor(df4$group, levels = c('安慰剂BID','2.5mgBID','5mgBID'))

lme_model4 <-
  lmerTest::lmer(score ~ time + group + W1D1 + drughis + (1|`PatientID`),
                 data = df4
  )

summary(lme_model4)
anova(lme_model4)

plot(ggemmeans(lme_model4, terms = c("time", "group"),
               condition = c(diagnose = "severe"))) +
  ggplot2::ggtitle("GLMER Effect plot")

