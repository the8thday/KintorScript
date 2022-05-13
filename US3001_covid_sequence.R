## ---------------------------
##
## Script name: title here
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2022-03-22
##
## Copyright (c) cliu, 2022
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: GT0918-US-3001测序SARS-COV-2的后续结果
##
##
## ---------------------------


library(tidyverse)


lineage1 <- read_csv('~/OneDrive/kintor/Daily_Work/US30001/GT0918-US-3001- Lineage Table-09Mar2022.csv') %>%
  janitor::remove_empty()
# lineage2 包含lineage1
lineage2 <- read_csv('~/OneDrive/kintor/Daily_Work/US30001/GT0918-US-3001- Lineage-table updated 15MAR2022.csv') %>%
  janitor::remove_empty()

variant <- read_csv('~/OneDrive/kintor/Daily_Work/US30001/GT0918-US-3001- Variant Table-09Mar2022.csv')

# 833例 实验室检测数据
blindbio <- read_csv('~/OneDrive/kintor/Daily_Work/US30001/GT0918US3001_BlindBio_14MAR2022.csv') %>%
  janitor::clean_names() %>%
  janitor::remove_empty()
# 828例 病毒载量数据
viralLoad <- read_csv('~/OneDrive/kintor/Daily_Work/US30001/GT0918US3001_ViralLoad_14MAR2022.csv') %>%
  janitor::clean_names() %>%
  janitor::remove_empty()


# lineage2 summary
skimr::skim(lineage2) %>% as_tibble()

# 统计病毒株系分布
lineage2 %>% janitor::tabyl(Lineage) %>%
  arrange(desc(n)) %>%
  filter(!Lineage %in% c('Not evaluable, will be repeated in the next batch', 'Not evaluable')) %>%
  ggpubr::ggbarplot(x = 'Lineage', y = 'n', label = T, fill = "steelblue", color = "steelblue") +
  theme(axis.text.x = element_text(angle = 90))
# 将AY的拼在一起的病毒系分布
lineage2 %>%
  mutate(Lineage2 = ifelse(str_detect(Lineage, '^AY'), 'Delta', Lineage)) %>%
  mutate(Lineage2 = ifelse(str_detect(Lineage, '^BA'), 'Omicron', Lineage2)) %>%
  mutate(Lineage2 = ifelse(Lineage2=='B.1.1.7', 'Alpha', Lineage2)) %>%
  mutate(Lineage2 = ifelse(Lineage2=='P.1', 'Gamma', Lineage2)) %>%
  mutate(Lineage2 = ifelse(Lineage2=='B.1.617.2', 'Delta', Lineage2)) %>%
  janitor::tabyl(Lineage2) %>%
  arrange(desc(n)) %>%
  filter(!Lineage2 %in% c('Not evaluable, will be repeated in the next batch', 'Not evaluable')) %>%
  ggpubr::ggbarplot(x = 'Lineage2', y = 'n', label = T, fill = "steelblue", color = "steelblue") +
  theme(axis.text.x = element_text(angle = 90))


# 各样本突变种类分布, 还只有150例的样本
variant %>% group_by(uniqueSampleID) %>% summarise(var_num = n()) %>%
  arrange(desc(var_num))


















