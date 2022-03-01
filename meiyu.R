## ---------------------------
##
## Script name: meiyu.R
##
## Purpose of script: pivot PK Data
##
## Author: LiuCong
##
## Date Created: 2022-01-19
##
## Copyright (c) cliu, 2022
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 美博所需之数据整理。包含两部分，第一部分为上海观合PK数据；第二部分为二次编码数据
##
##
## ---------------------------
# read files

library(tidyverse)
# library(tidyxl)
# library(magrittr)
library(ggprism)

#
# 二次编码 --------------------------------------------------------------------

pivot_data <- function(datapath, outpath, ...){
  stopifnot(file.exists(datapath))
  if(!dir.exists(outpath)){print('outpath doesnot exist!')}
  df <- readxl::read_excel(datapath,
                           sheet = 'Sheet1')

  doses <- df %>% select(dose) %>% distinct() %>% pull(dose)
  if(length(doses)==1){
    print('Only one dose class!')
    res1 <- df %>% filter(dose == doses[1]) %>% replace_na(list('Day Nominal'=0, `Hour Nominal`=1000)) %>%
      select(Subject, `Day Nominal`, `Hour Nominal`, `Concentration (pg/mL)`) %>%
      pivot_wider(names_from = Subject, values_from = `Concentration (pg/mL)`,
                  names_prefix = paste0(doses[1],'_')) %>%
      unite(col = 'DayHour', c(`Day Nominal`, `Hour Nominal`), remove = F)
    res1 %>% write_delim(file = file.path(outpath,'res.txt'), delim = '\t')
    return(res1)
  }else{
    res1 <- df %>% filter(dose == doses[1]) %>% replace_na(list(`Day Nominal`=0, `Hour Nominal`=1000)) %>%
      select(Subject, `Day Nominal`, `Hour Nominal`, `Concentration (pg/mL)`) %>%
      pivot_wider(names_from = Subject, values_from = `Concentration (pg/mL)`,
                  names_prefix = paste0(doses[1],'_')) %>%
      unite(col = 'DayHour', c(`Day Nominal`, `Hour Nominal`), remove = F)
    for(i in 2:length(doses)){
      print(paste0("Pivote sample ",doses[i]))
      res <- df %>% filter(dose == doses[i]) %>% replace_na(list(`Day Nominal`=0, `Hour Nominal`=1000)) %>%
        select(Subject, `Day Nominal`, `Hour Nominal`, `Concentration (pg/mL)`) %>%
        pivot_wider(names_from = Subject, values_from = `Concentration (pg/mL)`,
                    names_prefix = paste0(doses[i],'_')) %>%
        unite(col = 'DayHour', c(`Day Nominal`, `Hour Nominal`), remove = TRUE)
      # print(res)
      res1 <- res1 %>% full_join(res, by = 'DayHour')
    }
  }
  res1 %>% write_delim(file = file.path(outpath,'res.txt'), delim = '\t')
}


pivot_data2 <- function(datapath, outpath,
                        Subject = "Subject",
                        Day = "Day Nominal",
                        Hour = "Hour Nominal",
                        con = "Concentration (pg/mL)",
                        dose = "dose",
                        ...) {
  stopifnot(file.exists(datapath))
  if (!dir.exists(outpath)) {
    print("outpath doesnot exist! Try to create it!")
    dir.create(outpath, recursive = TRUE)
  }
  df <- readxl::read_excel(datapath,
    sheet = "Sheet1"
  )

  doses <- df %>%
    select(all_of(dose)) %>%
    distinct() %>%
    pull(all_of(dose))
  if (length(doses) == 1) {
    print("Only one dose class!")
    res1 <- df %>%
      filter(.data[[dose]] == doses[1]) %>%
      dplyr::mutate(across(all_of(Day), ~ replace_na(.x, 0))) %>%
      dplyr::mutate(across(all_of(Hour), ~ replace_na(.x, 1000))) %>%
      select(all_of(c(Subject, Day, Hour, con))) %>%
      pivot_wider(
        names_from = all_of(Subject), values_from = all_of(con),
        names_prefix = paste0(doses[1], "_")
      ) %>%
      unite(col = "DayHour", all_of(c(Day, Hour)), remove = F)
    res1 %>% write_delim(file = file.path(outpath, "res.txt"), delim = "\t")
    return(res1)
  } else {
    res1 <- df %>%
      filter(.data[[dose]] == doses[1]) %>%
      dplyr::mutate(across(all_of(Day), ~ replace_na(.x, 0))) %>%
      dplyr::mutate(across(all_of(Hour), ~ replace_na(.x, 1000))) %>%
      select(all_of(c(Subject, Day, Hour, con))) %>%
      pivot_wider(
        names_from = all_of(Subject), values_from = all_of(con),
        names_prefix = paste0(doses[1], "_")
      ) %>%
      unite(col = "DayHour", all_of(c(Day, Hour)), remove = F)
    for (i in 2:length(doses)) {
      print(paste0("Pivote sample ", doses[i]))
      res <- df %>%
        filter(.data[[dose]] == doses[i]) %>%
        dplyr::mutate(across(all_of(Day), ~ replace_na(.x, 0))) %>%
        dplyr::mutate(across(all_of(Hour), ~ replace_na(.x, 1000))) %>%
        select(all_of(c(Subject, Day, Hour, con))) %>%
        pivot_wider(
          names_from = all_of(Subject), values_from = all_of(con),
          names_prefix = paste0(doses[i], "_")
        ) %>%
        unite(col = "DayHour", all_of(c(Day, Hour)), remove = TRUE)
      # print(res)
      res1 <- res1 %>% full_join(res, by = "DayHour")
    }
  }
  res1 %>%
    arrange(.data[[Day]]) %>%
    write_delim(file = file.path(outpath, "res2.txt"), delim = "\t")
}

# 修改以下需要修改的地方
datapath <- '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/GT20029-DATA-20211117-二次编码.xlsx'
pivot_data(datapath, '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/')
pivot_data2(datapath,
            outpath = '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/',
            Subject = 'Subject',
            Day = 'Day Nominal',
            Hour = 'Hour Nominal',
            con = 'Concentration (pg/mL)',
            dose = 'dose')


# PK 检测 -------------------------------------------------------------------

extract_PK <- function(filepath, outpath,
                       subject = '受试者编号',
                       con = '检测浓度',
                       time = '理论采集时间点',
                       date = '研究日',
                       dose = '剂量组',
                       ...) {
  stopifnot(file.exists(filepath))
  if(!dir.exists(outpath)){
    print('outpath doesnot exist! Try to create it!')
    dir.create(outpath, recursive = TRUE)
  }

  df <- readxl::read_excel(filepath)
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

  df2 <- df %>% select(`研究日`, `理论采集时间点`, `受试者编号`, `剂量组`, `检测浓度`) %>%
    unite('key', c(`受试者编号`, `剂量组`)) %>%
    pivot_wider(names_from = key, values_from = `检测浓度`)
  df2 %>% write_delim(file.path(outpath,'pivot_data.txt'), delim = '\t')
}

extract_PK2 <- function(filepath, outpath,
                       subject = '受试者编号',
                       con = '检测浓度',
                       time = '理论采集时间点',
                       date = '研究日',
                       dose = '剂量组',
                       ...) {
  stopifnot(file.exists(filepath))
  if(!dir.exists(outpath)){
    print('outpath doesnot exist! Try to create it!')
    dir.create(outpath, recursive = TRUE)
  }

  df <- readxl::read_excel(filepath)
  print(table(df[[subject]]))

  df <- df %>%
    mutate(across(all_of(con), ~ifelse(.x=='BQL', 0, .x))) %>%
    mutate(across(all_of(con), ~as.numeric(.x))) %>%
    mutate(
      across(all_of(time), ~case_when(
        .data[[date]] == 'C0D2'~ '24h',
        .data[[date]] == 'C0D3'~ '48h',
        .data[[date]] == 'C0D4'~ '72h',
        TRUE ~ .x
      ))
    ) %>%
    mutate(across(all_of(time), ~ifelse(.x=='给药前', 0, .x))) %>%
    mutate(across(all_of(time), ~str_extract(.x, '\\d+(\\.\\d)?'))) %>%
    arrange(.data[[dose]])

  df2 <- df %>% select(all_of(c(date, time, subject, dose, con))) %>%
    unite('key', all_of(c(subject, dose))) %>%
    pivot_wider(names_from = key, values_from = all_of(con))
  df2 %>% write_delim(file.path(outpath,'pivot_data.txt'), delim = '\t')
}

inpath <- '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/GT0486-CN-1001_ PK检测_上海观合医药科技有限公司_非机密_20210526_01.xls'
# dd <- extract_PK(inpath, '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/')

dd <- extract_PK2(inpath, '/Users/congliu/OneDrive/kintor/Daily_Work/meiyu/',
           subject = '受试者编号',
           con = '检测浓度',
           time = '理论采集时间点',
           date = '研究日',
           dose = '剂量组')

# 求取均值因为分类数目不等
dd %<>% rowwise() %>%
  mutate(avg_40mg = mean(c_across(ends_with('40mg'))))


dd %>%
  select(-ends_with('40mg')) %>%
  filter(`研究日` %in% c('C0D1','C1D28')) %>%
  mutate(`理论采集时间点`=as.numeric(`理论采集时间点`)) %>%
  pivot_longer(-c(`研究日`, `理论采集时间点`)) %>%
  unite('key', c(name, `研究日`), remove = F) %>%
  ggplot(., aes(x = `理论采集时间点`, y = value, color=key, group = key)) +
  geom_point(aes(shape = key), size=4) +
  geom_line(size=1) +
  ggsci::scale_color_aaas() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top")
        ) + xlab('Time(hr)') + ylab('Plasma Conc(ng/ml)')

dd$研究日 <- factor(dd$研究日, levels = unique(dd$研究日))

dd %>% filter(!str_detect(`研究日`, 'C0')) %>%
  select(c(`研究日`, `2001_10mg`, `1001_20mg`)) %>%
  pivot_longer(-`研究日`) %>%
  ggplot(., aes(x = `研究日`, y = value, color=name, group = name)) +
  geom_point(aes(shape = name), size=4) +
  geom_line(size=1) +
  ggsci::scale_color_aaas() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top")
  ) + xlab('Days') + ylab('Plasma Conc(ng/ml)')


