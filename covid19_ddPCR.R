## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2022-02-14
##
## Copyright (c) cliu, 2022
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 一批次没有揭盲的新冠患者的ddPCR检测结果。主要是一些描述性的统计分析。新冠病毒载量的变化数据。
##
##
## ---------------------------


library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)


ddPCR <- read_csv('/Users/congliu/OneDrive/kintor/Daily_Work/covid19_ddpcr/ddPCR Viral load report_Feb 11.csv') %>%
  janitor::clean_names() %>% janitor::remove_empty()

# 共745例样本； TestName 共有9种，不过排除Repeat ID的共有6种检测；Visit Name为SCREEN/D3/D7/D14/D28排除EW/UNSCHED
janitor::tabyl(ddPCR$test_name)


test_name2 <- 'COVID ddPCR NP final result'

ddPCR %>% filter(test_name == test_name2) %>%
  filter(!patient_id %in% abnormal$patient_id) %>%
  janitor::tabyl(visit_name, test_result) %>%
  janitor::adorn_totals(where = 'col') %>%
  knitr::kable()

test_res <- ddPCR %>% filter(test_name == test_name2) %>%
  filter(!visit_name %in% c('EW', 'UNSCHED')) %>%
  # filter(!`Test Result` == 'Detected but below assay LoQ') %>%
  select(`patient_id`, visit_name, `test_name`, test_result) %>%
  pivot_wider(names_from = visit_name, values_from = `test_result`)

test_res

# 以热图的形式展示不同TestName的不同TestResult
M <- test_res %>% select(-test_name) %>%
  filter(!patient_id %in% abnormal$patient_id) %>%
  column_to_rownames('patient_id')

colors <- structure(c("red", "blue",'green', 'black'),
                    names = c('Positive', 'Negative', 'Detected but below assay LoQ','NOT REPORTED'))

Heatmap(M,
        name = test_name2,
        na_col = 'white',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        col = colors,
        column_names_rot = 90,
        show_column_names = TRUE,
        show_row_names = FALSE,
        width = unit(2.5, 'cm'),
        row_names_gp = gpar(fontsize=8)
        )
# 挑选异常样本
abnormal <- test_res %>% filter((SCREEN == 'Negative' & (D3 != 'Negative' | D7 != 'Negative' | D14 != 'Negative' | D28 != 'Negative')) |
                      (D3 == 'Negative' & (D7 != 'Negative' | D14 != 'Negative' | D28 != 'Negative')) |
                      (D7 == 'Negative' & (D14 != 'Negative' | D28 != 'Negative')) |
                      (D14 == 'Negative' & D28 != 'Negative')
                    ) %>% distinct() %>%
  relocate(D7, .before = D14)
M <- abnormal %>% select(-`test_name`) %>%
  column_to_rownames('patient_id')


## 针对连续变量的探索
# 首先看一下此结局变量的分布
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

test_name1 <- 'N2 copies/mL (Np)'

test_res3 <- ddPCR %>% filter(test_name == test_name1) %>%
  filter(!visit_name %in% c('EW', 'UNSCHED')) %>%
  select(patient_id, visit_name, test_name, test_result, `investigator_name`) %>%
  mutate(test_result = as.numeric(test_result),
         patient_id = as.character(patient_id)
  )

ggpubr::gghistogram(test_res3, x = 'test_result')
ggpubr::ggqqplot(test_res3, x = 'test_result') # 数据严重偏态

test_res3 <- test_res3 %>%
  mutate(test_result  = log(test_result +1))
test_res3$visit_name <- factor(test_res3$visit_name, levels = c('SCREEN', 'D3', 'D7', 'D14', 'D28'))

test_res3 %>% filter(is.na(test_result)) %>% janitor::tabyl(visit_name)

ggplot(data = test_res3, aes(x = `investigator_name`, y = `test_result`, color = visit_name)) +
  geom_boxplot(outlier.size = 1) +
  ggsci::scale_color_aaas() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

# summary analysis
test_res2 <- ddPCR %>% filter(test_name == test_name1) %>%
  filter(!visit_name %in% c('EW', 'UNSCHED')) %>%
  select(patient_id, visit_name, patient_id, test_result) %>%
  pivot_wider(names_from = visit_name, values_from = test_result) %>%
  mutate(across(starts_with('D'), as.numeric)) %>%
  mutate(SCREEN = as.numeric(SCREEN),
         patient_id = as.character(patient_id)
         )

test_res2
visdat::vis_miss(test_res2)
gtsummary::tbl_summary(test_res2 %>% select(-patient_id),
                       type = all_continuous() ~ "continuous2",
                       statistic = all_continuous() ~ c("{N_nonmiss}",
                                                        "{median} ({p25}, {p75})",
                                                        "{min}, {max}")
                       )








