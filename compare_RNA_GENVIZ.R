library(tidyverse)


# compare DEGs ------------------------------------------------------------

path <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Mouse/Report/Result/08_DE'

all_fpkm <- read_delim(
  '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Mouse/Report/Result/06_GeneExpression/all.fpkm_anno.xls',
  delim = '\t'
  )


R1R2 <- read_delim(file.path(path, 'Group_R-1-VS-R-2_DE_anno.xls'),
                   delim = '\t')

R1R2_sig <- read_delim(file.path(path, 'Group_R-1-VS-R-2_DE_significant_anno.xls'),
                       delim = '\t')
R2R3_sig <- read_delim(file.path(path, 'Group_R-2-VS-R-3_DE_significant_anno.xls'),
                       delim = '\t')
R2R4_sig <- read_delim(file.path(path, 'Group_R-2-VS-R-4_DE_significant_anno.xls'),
                       delim = '\t')


go_inter <- function(x){
  x %>%
    select(geneID, pval, FDR, logFC, Regulation, GeneSymbol)
}


R1R2_R2R3 <- go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R3_sig), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  select(sort(names(.)))

R1R2_R2R3 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/R1R2_R2R3.txt',
                          delim = '\t')


go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R4_sig), by = 'geneID', suffix = c('_R1R2', '_R2R4')) %>%
  select(sort(names(.))) %>%
  write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/R1R2_R2R4.txt', delim = '\t')


foo <- go_inter(R2R4_sig) %>%
  rename_with(.fn = ~paste0(.x, '_R2R4')) %>%
  rename(geneID=geneID_R2R4)

go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R3_sig), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  full_join(foo, by = 'geneID') %>%
  select(sort(names(.)))


# GO ----------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)





