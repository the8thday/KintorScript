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
library(lme4)
# library(WRS2) # to keep the outliers in the data and perform robust ANOVA test using the WRS2 package


hormones <- c('ESR', 'DHT', 'testos', 'Estradiol', 'SHBG')
sex_list <- c('males', 'female')
jisu_sheet <- 'ESR'
sex <- 'males'

mt_ESR <- readxl::read_excel(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_treated.xlsx'),
                             sheet = jisu_sheet
                             ) %>%
  select(-BMI) %>%
  drop_na()
mn_ESR <- readxl::read_excel(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_untreated.xlsx'),
                             sheet = jisu_sheet
                             ) %>%
  select(-BMI) %>%
  drop_na()

# 单变量重复资料 -----------------------------------------------------------------

long_data <- mn_ESR %>%
  rename('PatientID'='Patient ID') %>%
  pivot_longer(cols = -c(`PatientID`, AGE),
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
(get_anova_table(res.aov) %>%
  write_delim(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_{jisu_sheet}_onewayANOVA.txt'),
              delim = '\t'
              ))

# use aov & post-hoc analysis
aov.res <- aov(score ~ time + Error(`PatientID`/time),
               data = long_data)
report::report(aov.res)
summary(aov.res)

pwc <- long_data %>%
  pairwise_t_test(
    score ~ time,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )
(pwc %>%
  write_delim(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_{jisu_sheet}_onewayANOVA_post.txt'),
              delim = '\t'
  ))

# friedman test
(res.fried <- long_data %>% friedman_test(score ~ time | PatientID))
res.fried %>%
  write_delim(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_{jisu_sheet}_friedman.txt'),
              delim = '\t'
  )
pwc2 <- long_data %>%
  wilcox_test(score ~ time, paired = TRUE,
              p.adjust.method = "bonferroni")
(pwc2 %>%
    write_delim(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_{jisu_sheet}_friedman_post.txt'),
                delim = '\t'
    ))

# 单因素协方差, ANCOVA
aov_xie <- aov(score ~ AGE + time, data = long_data) # 协变量写在因子前面
summary(aov_xie)
car::Anova(aov_xie, type = 'III')
postHocs <- multcomp::glht(aov_xie, linfct = multcomp::mcp(time = "Tukey"))
summary(postHocs)
sm::sm.ancova(x = long_data$AGE, y = long_data$score, group = long_data$time)

# plot
long_data %>%  ggplot(aes(x = time, y = score, color = time, fill = time)) +
  geom_boxplot(na.rm = TRUE,
               # outlier.shape = NA
               ) +
  ggtitle(glue::glue('{sex} {jisu_sheet}')) +
  theme_prism(base_size = 14) +
  scale_colour_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  geom_signif(
    y_position = c(130, 130, 130), xmin = c(0.8, 0.8, 2.8), xmax = c(1.2, 3.2, 3.2),
    annotation = c("*", "NS", "NS"), tip_length = 0
  )
  geom_signif(
    comparisons = list(c("Day 0", "Day 1"), c("Day 0", "Day 7"), c('Day 1', 'Day 7')),
    map_signif_level = TRUE,
    textsize = 4,
    test = "wilcox.test",
    test.args = list(paired = TRUE),
    vjust = 0.2
  )

ggwithinstats(data = long_data,
              x = time,
              y = score,
              type ='nonparametric'
               )

ggboxplot(
  data = long_data,
  x = 'time',
  y = 'score',
  fill = 'time',
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
) + stat_compare_means(comparisons =list(c("Day 0", "Day 1"), c("Day 0", "Day 7"), c('Day 1', 'Day 7')),
                       method = 'wilcox.test',
                       paired = TRUE,
                       aes(label = "p.signif")
                       )+
  stat_compare_means(label.y =101)


#
# two.way 重复资料 ------------------------------------------------------------

# plot bar, time为组内变量、class为组间变量
df_bar <- mt_ESR %>%
  mutate(class = 'Treated') %>%
  bind_rows(
    mn_ESR %>% mutate(class = 'Untreated')
  ) %>%
  pivot_longer(cols = -c(`Patient ID`, class, AGE), names_to = 'time', values_to = 'score') %>%
  convert_as_factor(class, time) %>%
  drop_na()

df_bar %>%
  group_by(class, time) %>%
  shapiro_test(score)

two.way <- anova_test(score ~ time * class + Error(`Patient ID`/time),
                      data = df_bar
)
(get_anova_table(two.way) %>%
  write_delim(glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/{sex}_{jisu_sheet}_twowayANOVA.txt'),
              delim = '\t'
  ))
# post-hoc
pwc <- df_bar %>%
  group_by(class) %>%
  pairwise_t_test(
    score ~ time,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc


res.aov <- anova_test(
  data = df_bar,
  dv = score, wid = `Patient ID`,
  within = c(time),
  between = class,
  covariate = AGE
)
get_anova_table(res.aov)

fit <- aov(score ~ class * time + Error(`Patient ID`/time),
           data = df_bar)
# report::report(fit)
summary(fit)
anova_summary(fit)

# 协方差分析ANCOVA
aov_xie2 <- aov(score ~ AGE + time, data = df_bar) # 协变量写在因子前面
summary(aov_xie2)
car::Anova(aov_xie2, type = 'III')

# 广义估计方程 GEE
geefit <- geeglm(score ~ time + class,
               id=`Patient ID`,
               corstr='exchangeable',
               family="gaussian",
               data=df_bar,
               std.err = 'san.se')
summary(geefit)
broom::tidy(geefit)

# 线性混合模型
model_lmer <- lme4::lmer(score ~ time + class + (1|`Patient ID`) + (1+time|`Patient ID`),
                         data = df_bar
)
performance::check_model(model_lmer)
summary(model_lmer)


three_p <- df_bar %>% group_by(time) %>%
  wilcox_test(score ~ class) %>% pull(p)

print(c(glue::glue('p_0 P value: {three_p[1]}'),
        glue::glue('p_1 P value: {three_p[2]}'),
        glue::glue('p_7 P value: {three_p[3]}')
        ))

df_bar %>%
  group_by(time, class) %>%
  summarise(count = n(),
            mean  = mean(score, na.rm = TRUE),
            median = median(score, na.rm = TRUE),
            sd = sd(score, na.rm = TRUE),
            N = length(score),
            se = sd/sqrt(N)
            ) %>%
  ggplot(data = ., mapping = aes(x = time, y = mean, fill = class)) +
  geom_col(position = "dodge", width = .5) +
  geom_errorbar(mapping = aes(
                              ymin = mean-se,
                              ymax = mean + se,
                              width = 0.2),
                position = position_dodge(0.5)
                ) +
  geom_signif(
    y_position = c(80, 80, 80), xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2),
    annotation = c("NS", "NS", "NS"), tip_length = 0
  ) +
  theme_prism(base_size = 14) +
  scale_fill_prism(palette = "floral") +
  ggtitle(glue::glue('{jisu_sheet} Treated VS Untreated'))


ggpubr::ggbarplot(data = df_bar,
                  x = 'time',
                  y = 'score',
                  fill = 'class',
                  add = 'mean_se',
                  palette = 'jco',
                  position = position_dodge()
                  ) +
  stat_compare_means(comparisons =list(c("Day 0", "Day 1"), c("Day 0", "Day 7"), c('Day 1', 'Day 7')),
                                         method = 'wilcox.test',
                                         paired = TRUE,
                                         aes(label = "p.signif")
                  )








