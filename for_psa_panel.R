library(tidyverse)
library(maftools)



cfdna_snv <- readxl::read_excel('D:kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_202010120.xlsx',
                                sheet = 'cfDNA_SNV')
# maf <- read.maf()
cfdna_cnv <- readxl::read_excel('D:kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_202010120.xlsx',
                                sheet = 'cfDNA_CNV')
cfdna_qc <- readxl::read_excel('D:kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_202010120.xlsx',
                               sheet = 'cfDNA_QC_Summary')
cfrna_splicing <- readxl::read_excel('D:kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_202010120.xlsx',
                                     sheet = 'cfRNA_Splicing')

metadata <- readxl::read_excel('D:kintor/2020-1014-general_dates.xlsx')
d <- metadata %>% 
  mutate(time = `EOT`-`Date of 1st Dosing`) %>% 
  mutate(stat = if_else(`Treatment discontunation`%in% 
                          c('Bone or Radiographic Disease Progression')|
                          `if others, specify`%in% c('Rise in PSA',
                                                     'PI decision due to rising PSA'
                          ), 
                        2, 1)) %>% 
  select(c('pt. ID', 'dosing', 'Prior Treat', 'time', 'stat')) %>% 
  rename(drug = `Prior Treat`)
d$time %<>% as.vector

df <- cfdna_snv %>% left_join(d, by = c('SubjectID'='pt. ID'))

df %>% mutate(time_d = ifelse(StudyVisit=='C1D1', time-0, ifelse(
  StudyVisit=='C4D1', time-84, ifelse(
    StudyVisit=='C7D1', time-84*2, ifelse(
      StudyVisit=='C10D1', time-84*3, ifelse(
        StudyVisit=='C13D1', time-84*4, ifelse(
          StudyVisit=='C16D1', time-84*5, ifelse(
            StudyVisit=='End of Treatment', 0, NA
          )
        )
      )
    )
  )
))) %>% select(SubjectID, ExternalID, time, time_d) %>% arrange(time_d)


# AR ----------------------------------------------------------------------
# 对于ARcnv而言
df2 <- cfdna_cnv %>% left_join(d, by = c('SubjectID'='pt. ID')) #对于发生splicing的样本而言
# df2 <- df %>% left_join(cfdna_cnv, by = c('SubjectID'='SubjectID'))

df2 %>% filter(SubjectID %in% long) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(SubjectID) %>% count()

df2 %>% filter(SubjectID %in% short) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(`Variant_Classification`) %>% 
  # count()
  select(ExternalID) %>% arrange(ExternalID)

# 对于AR splicing而言
df3 <- cfrna_splicing %>% left_join(d, by = c('SubjectID'='pt. ID'))

df3 %>% filter(SubjectID %in% long) %>% 
  group_by(SubjectID) %>% count()

df3 %>% filter(SubjectID %in% long) %>% group_by(NameAST) %>% 
  select(ExternalID) %>% arrange(ExternalID)


# AR-V7
df3 %>% filter(NameAST == 'AR-V7') %>% filter(SubjectID %in% short) %>% 
  distinct(SubjectID)

df3 %>% filter(NameAST =='AR-V7') %>% 
  filter(SupportReads >= 3) %>% 
  select(ExternalID, StudyVisit, time) %>% arrange(desc(time))



# tp53 --------------------------------------------------------------------
# 就tp53的分布而言，呈现出药效短的人其TP53突变较高, 进一步分析TP53的突变类型



# summary analysis --------------------------------------------------------

cfdna_snv %>% group_by(SubjectID) %>% count()
# 查看每个样本所对应的时期的总数，采样时间的分布
foo <- cfdna_qc %>% group_by(SubjectID) %>% count() %>% arrange(desc(n))
ggplot(foo) +
  geom_bar(
    aes(x = n)
  )
cfdna_qc %>% group_by(SubjectID) %>% 
  filter(n()>2)
# 进一步查看每个时期的分析
cfdna_qc




# 以30号和1号为例
# 设时间长短的两类，分别按照变异数目、变异类别、变异频率、CNV、AR
# 
bar <- df %>% select(c(SubjectID, ExternalID, StudyVisit, Variant_Value, 
                       Hugo_Symbol, HGVSc, HGVSp, time, stat)) %>% 
  arrange(desc(time)) %>% 
  unite('variant', c(Hugo_Symbol, HGVSp), sep = ':', remove = T)

# 41个样本按照PFS时间排序
t <- bar %>% select(SubjectID, time) %>% distinct() %>% drop_na()
long <- t$SubjectID[1:11]
short <- t$SubjectID[12:32]
df %>% filter(SubjectID %in% long) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(SubjectID) %>% count()
# 21个short样本中TP53突变的分布情况
df %>% filter(SubjectID %in% long) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(HGVSp_Short) %>% 
  # count()
  select(ExternalID, Clinvar, CLIN_SIG) %>% arrange(ExternalID)












