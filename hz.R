## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: cliu
##
## Date Created: 2021-06-25
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


inpath <- '/Users/congliu/OneDrive/kintor/相反数一对一.xlsx'
outpath <- '/Users/congliu/OneDrive/kintor/'
cohortname <- 'test'

hz <- readxl::read_excel(inpath) %>% 
  select(-remark) %>% 
  arrange(item)

find_index <- function(aa_v){
  if(length(aa_v)==1){return(0)}
  stop <- FALSE
  for(i in 1:(length(aa_v)-1)){
    for(j in (i+1):length(aa_v)){
      if(aa_v[i]+aa_v[j]==0){
        f1 <- i
        f2 <- j
        stop <- TRUE
        break
      }
    }
    if(stop){break}
  }
  c(f1, f2)
}


nest.index <- hz %>% nest_by(item) %>% 
  mutate(index = list(find_index(data$amount)))

key.index <- hz %>% group_by(item) %>% 
  summarise(index=find_index(.data$amount)) %>% 
  filter(index != 0) %>% 
  ungroup() %>% 
  unite('key',item:index, remove = F) %>% 
  select(-item)


final_result <- hz %>% group_by(item) %>% 
  mutate(ob = 1:n()) %>% 
  ungroup() %>% 
  unite('key', c(item,ob), remove = F) %>% 
  left_join(key.index, by = c('key')) %>% 
  mutate(remark=if_else(
    is.na(index), 'NA', 'offset'
  )) %>% 
  select(-c(ob,index,key))

write_delim(final_result, 
            file.path(
              outpath,
              paste0("remark_offset", cohortname, ".txt")
            ), delim = '\t')





