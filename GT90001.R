
library(tidyverse)


(blood <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/Daily_Work/GT90001/blood chemistry_2.xlsx'
  ) %>%
    filter(AnalyteName=='Triglyceride (TG)') %>%
    # filter(AnalyteName=='Total cholesterol\r\n(CHOL)') %>%
  select(Subject, Period, AnalyteValue) %>%
  filter(AnalyteValue != 'ND') %>%
  filter(!str_detect(Period, 'Unscheduled')) %>%
  mutate(AnalyteValue = as.numeric(AnalyteValue),
         Period = gsub(' \\(1\\)', '', Period),
         Period = gsub('\r\n', ' ', Period)
         ) %>% rstatix::convert_as_factor(Period) %>% filter(Period != 'Post-treatment follow up')
    # filter(!Period %in% c('Cycle 14 Day 01', 'Cycle 34 Day 01', 'Cycle 36 Day 01'))
  # filter(!(Subject == 'G01013' & Period == 'Cycle 03 Day 01'))
  )
blood$Period <- factor(blood$Period, levels = levels(fct_relevel(blood$Period, '筛选中')))

table(blood$Subject, blood$Period) # 21例患者的55个时期的数据



# plot 所有样本的均值
ggplot(blood, mapping = aes(x = Period, y = AnalyteValue, group = 1)) +
  stat_summary(geom = 'line',
               fun = 'mean', color = 'red'
  ) +
  stat_summary(geom = 'point',
               fun = 'mean', color = 'blue'
  ) +
  stat_summary(geom = 'errorbar',
               fun.min = function(x){mean(x)-sd(x)},
               fun.max = function(x){mean(x)+sd(x)},
               width = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0, NA)
  ) +
  ggsci::scale_color_lancet() +
  ggprism::theme_prism() +
  theme(text = element_text(family='STHeiti'),
        axis.text.x = element_text(angle = 90, size = 6)
        )


ggplot(blood,
       aes(x = Period, y = AnalyteValue, color = Subject)) +
  geom_point() +
  geom_line(aes(group = Subject)) +
  # facet_grid(Subject ~ ., scales = 'free') +
  scale_y_continuous(expand = expansion(mult = c(0,0)),
                     limits = c(0, NA)
  ) +
  # ggsci::scale_color_lancet() +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(21)) +
  ggprism::theme_prism() +
  # theme_bw() +
  theme(text = element_text(family='STHeiti'),
        axis.text.x = element_text(angle = 90, size = 6),
        legend.position="bottom"
        # axis.text.y = element_text(size = 4),
        # strip.text.y = element_text(angle = 0, size = 6),
        # legend.position = 'none'
  )




