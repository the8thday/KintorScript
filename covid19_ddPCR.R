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
## Notes:
##
##
## ---------------------------


library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)


ddPCR <- read_csv('/Users/congliu/OneDrive/kintor/Daily_Work/covid19_ddpcr/ddPCR Viral load report_Feb 11.csv')

# 共745例样本； TestName 共有9种，不过排除Repeat ID的共有6种检测；Visit Name为SCREEN/D3/D7/D14/D28排除EW/UNSCHED
table(ddPCR$`Test Name`)

test_name <- 'COVID ddPCR NP final result'
test_res <- ddPCR %>% filter(`Test Name` == test_name) %>%
  filter(!`Visit Name` %in% c('EW', 'UNSCHED')) %>%
  select(`Patient ID`, `Visit Name`, `Test Name`, `Test Result`) %>%
  pivot_wider(names_from = `Visit Name`, values_from = `Test Result`)

test_res

# 以热图的形式展示不同TestName的不同TestResult
M <- test_res %>% select(-`Test Name`) %>%
  column_to_rownames('Patient ID')

colors <- structure(c("red", "blue",'green', 'black'),
                    names = c('Positive', 'Negative', 'Detected but below assay LoQ','NOT REPORTED'))

Heatmap(M,
        name = test_name,
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
        width = unit(2.5, 'cm')
        # column_names_gp = gpar(fontsize=1)
        )

test_name <- 'N1 copies/mL (Np)'
test_res2 <- ddPCR %>% filter(`Test Name` == test_name) %>%
  filter(!`Visit Name` %in% c('EW', 'UNSCHED')) %>%
  select(`Patient ID`, `Visit Name`, `Test Name`, `Test Result`) %>%
  pivot_wider(names_from = `Visit Name`, values_from = `Test Result`) %>%
  mutate(across(starts_with('D'), as.numeric)) %>%
  mutate(SCREEN = as.numeric(SCREEN),
         `Patient ID` = as.character(`Patient ID`)
         )

test_res2

gtsummary::tbl_summary(test_res2 %>% select(-`Patient ID`),
                       type = all_continuous() ~ "continuous2",
                       statistic = all_continuous() ~ c("{N_nonmiss}",
                                                        "{median} ({p25}, {p75})",
                                                        "{min}, {max}")
                       )





