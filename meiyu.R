
# read files

library(tidyverse)
library(magrittr)
library(ggprism)


df <- readxl::read_excel(
  'D:/GT0486-CN-1001_ PK检测_上海观合医药科技有限公司_非机密_20210526_01.xls')

print(table(df$受试者编号))


df <- df %>% 
  mutate(`检测浓度`=ifelse(`检测浓度`=='BQL',0,`检测浓度`)) %>% 
  mutate(`检测浓度`=as.numeric(`检测浓度`)) %>% 
  mutate(
    `理论采集时间点`=case_when(
      `研究日` == 'C0D2'~ '24h',
      `研究日` == 'C0D3'~ '48h',
      `研究日` == 'C0D4'~ '72h',
      TRUE ~ `理论采集时间点`
    )
  ) %>% 
  mutate(`理论采集时间点`=ifelse(`理论采集时间点`=='给药前', 
                          0, `理论采集时间点`)) %>% 
  mutate(`理论采集时间点`=str_extract(`理论采集时间点`,'\\d+(\\.\\d)?')) %>%
  arrange(`剂量组`)

dd <- df %>% unite('key', c('受试者编号','剂量组'), remove=F) %>% 
  select(key, `研究日`, `检测浓度`, `理论采集时间点`) %>% 
  pivot_wider(names_from = key, values_from = `检测浓度`)
  
dd %>% write_delim('', delim = '\t')

# 求取均值因为分类数目不等
dd %<>% rowwise() %>% 
  mutate(avg_40 = mean(c_across(ends_with('40mg'))))


dd %>% 
  select(-ends_with('40mg')) %>% 
  filter(`研究日` %in% c('C0D1','C1D28')) %>% 
  mutate(`理论采集时间点`=as.numeric(`理论采集时间点`)) %>% 
  pivot_longer(-c(`研究日`, `理论采集时间点`)) %>% 
  unite('key', c(name, `研究日`), remove = F) %>% 
  ggplot(., aes(x = `理论采集时间点`, y = value, color=key, group = key)) +
  geom_point(size=4, shape = 17) +
  geom_line(size=1) +
  ggsci::scale_color_aaas() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top")
        )


df %>% filter(`理论采集时间点`==0) %>% 
  filter(`受试者编号` %in% c('2001','1001')) %>% 
  fct_reorder()
  ggplot(., aes(x = `研究日`, y = `检测浓度`, color=factor(`受试者编号`), 
                group = factor(`受试者编号`))) +
  geom_point(size=4, shape = 17) +
  geom_line(size=1) +
  ggsci::scale_color_aaas() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top")
  )





