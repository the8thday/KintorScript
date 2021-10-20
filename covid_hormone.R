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
## Notes:
##
##
## ---------------------------


library(tidyverse)
library(ggprism)
library(rstatix)
library(ggsignif)


jisu_sheet <- 'ESR'

mt_ESR <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/males_treated.xlsx',
                             sheet = jisu_sheet
                             )
mn_ESR <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/males_un.xlsx',
                             sheet = jisu_sheet
                             ) %>%
  drop_na()

mt_ESR %>%
  pivot_longer(cols = -`Patient ID`) %>%
  # anova_test(value ~ name)
  # kruskal_test(value ~ name)
  ggplot(aes(x = name, y = value, color = name, fill = name)) +
  geom_boxplot() +
  theme_prism(base_size = 14) +
  scale_colour_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  geom_signif(
    comparisons = list(c("Day 0", "Day 1"), c("Day 0", "Day 7"), c('Day 1', 'Day 7')),
    map_signif_level = TRUE, textsize = 6,
    test = "wilcox.test",
    vjust = 0.2
  )


df_bar <- mt_ESR %>%
  mutate(class = 'treat') %>%
  bind_rows(
    mn_ESR %>% mutate(class = 'untreated')
  ) %>%
  pivot_longer(cols = -c(`Patient ID`, class))

p_0 <- df_bar %>% filter(name == 'Day 0') %>%
  t_test(value ~ class) %>%
  pull(p)
p_1 <- df_bar %>% filter(name == 'Day 1') %>%
  t_test(value ~ class) %>%
  pull(p)
p_7 <- df_bar %>% filter(name == 'Day 7') %>%
  t_test(value ~ class) %>%
  pull(p)
print(c(glue::glue('p_0 P value: {p_0}'),
        glue::glue('p_1 P value: {p_1}'),
        glue::glue('p_7 P value: {p_7}')
        ))

df_bar %>%
  ggplot(data = ., mapping = aes(x = name, y = value, fill = class)) +
  geom_bar(stat = 'identity', position = "dodge", width = .5) +
  geom_signif(
    y_position = c(78, 65, 90), xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2),
    annotation = c("NS", "*", "**"), tip_length = 0
  ) +
  theme_prism(base_size = 14) +
  scale_fill_prism(palette = "floral")








