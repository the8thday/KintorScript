## ---------------------------
##
## Script name: covid_hormone.R
##
## Purpose of script: statistic covid patients hormone status
##
## Author: LiuCong
##
## Date Created: 2021-10-20
##
## Copyright (c) cliu, 2021
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 新冠患者在使用普克鲁胺后的一些激素水平在不同检测时间的变化
## 主要采用重复测量方差分析，
## 或将增加广义估计方程、线性混合模型分析、协方差分析
## ---------------------------


library(tidyverse)
library(ggprism)
library(rstatix)
library(ggsignif)
library(ggstatsplot)
library(geepack)


hormones <- c('ESR', 'DHT', 'testos', 'Estradiol', 'SHBG')
jisu_sheet <- 'testos'

mt_ESR <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/males_treated.xlsx',
                             sheet = jisu_sheet
                             ) %>% drop_na()
mn_ESR <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/males_untreated.xlsx',
                             sheet = jisu_sheet
                             ) %>% drop_na()

long_data <- mt_ESR %>%
  rename('PatientID'='Patient ID') %>%
  pivot_longer(cols = -`PatientID`,
               names_to = 'time',
               values_to = 'score'
               ) %>%
  drop_na() %>%
  convert_as_factor(PatientID, time)

# 检验是否符合anovar分析
long_data %>% group_by(time) %>% identify_outliers(score)
long_data %>% group_by(time) %>% shapiro_test(score)
long_data %>% levene_test(score ~ time)

aa <- long_data %>%
  anova_test(score ~ time + Error(`PatientID`/time)) # one way repeated measure ANOVA
res.aov <- anova_test(data = long_data, dv = score, wid = PatientID, within = time)
get_anova_table(res.aov)

# use aov & post-hoc analysis
aov.res <- aov(score ~ time + Error(`PatientID`/time), data = long_data)
report::report(aov.res)
summary(aov.res)
aov.res %>% tukey_hsd()

pwc <- long_data %>%
  pairwise_t_test(
    score ~ time, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

# friedman test
res.fried <- long_data %>% friedman_test(score ~ time |PatientID)
res.fried
pwc <- long_data %>%
  wilcox_test(score ~ time, paired = TRUE,
              p.adjust.method = "bonferroni")
pwc

# 协变量分析, 单因素协方差, ANCOVA
aov(score ~ age + time, data = long_data) # 协变量写在因子前面


long_data %>%  ggplot(aes(x = time, y = score, color = time, fill = time)) +
  geom_boxplot(na.rm = TRUE) +
  theme_prism(base_size = 14) +
  scale_colour_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  geom_signif(
    comparisons = list(c("Day 0", "Day 1"), c("Day 0", "Day 7"), c('Day 1', 'Day 7')),
    map_signif_level = TRUE, textsize = 6,
    test = "wilcox.test",
    # test.args = ''
    vjust = 0.2
  )

ggwithinstats(data = long_data,
              x = time,
              y = score,
              type ='nonparametric'
               )


# plot bar, name为组内变量、class为组间变量
df_bar <- mt_ESR %>%
  mutate(class = 'treat') %>%
  bind_rows(
    mn_ESR %>% mutate(class = 'untreated')
  ) %>%
  pivot_longer(cols = -c(`Patient ID`, class)) %>%
  drop_na()

df_bar %>%
  group_by(class, name) %>%
  shapiro_test(value)

two.way <- anova_test(value ~ name * class + Error(`Patient ID`/name),
                      data = df_bar
)
summary(two.way)

fit <- aov(value ~ name * class + Error(`Patient ID`/name),
           data = df_bar)
report::report(fit)
summary(fit)

# 协方差分析ANCOVA


# 广义估计方程 GEE
fit1 <- geeglm(value ~ name + class,
               id=`Patient ID`,
               corstr='exchangeable',
               family="gaussian",
               data=df_bar,
               std.err = 'san.se')
summary(fit1)

# 线性混合模型
model_lmer <- lme4::lmer(value ~ name * class + (1|cyl),
                         data = df_bar
)
performance::check_model(model_lmer)


three_p <- df_bar %>% group_by(name) %>%
  t_test(value ~ class) %>% pull(p)

print(c(glue::glue('p_0 P value: {three_p[1]}'),
        glue::glue('p_1 P value: {three_p[2]}'),
        glue::glue('p_7 P value: {three_p[3]}')
        ))

df_bar %>%
  group_by(name, class) %>%
  summarise(count = n(),
            mean  = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            N = length(value),
            se = sd/sqrt(N)
            ) %>%
  ggplot(data = ., mapping = aes(x = name, y = mean, fill = class)) +
  geom_col(position = "dodge", width = .5) +
  geom_errorbar(mapping = aes(x = name,
                              ymin = mean-se,
                              ymax = mean + se,
                              group = name,
                              width = 0.2),
                # position = position_dodge(width = 0.8)
                ) +
  geom_signif(
    y_position = c(78, 65, 90), xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2),
    annotation = c("NS", "*", "**"), tip_length = 0
  ) +
  theme_prism(base_size = 14) +
  scale_fill_prism(palette = "floral")


ggpubr::ggbarplot(data = df_bar,
                  x = 'name',
                  y = 'value',
                  fill = 'class',
                  add = 'mean_se',
                  palette = 'jco',
                  position = position_dodge()
                  )








