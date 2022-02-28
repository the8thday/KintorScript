# analysis for 1003

library(survival)
library(survminer)
library(tidyverse)
library(lubridate)
library(magrittr)

# meta1003 <- readxl::read_xlsx('D:/kintor/CN-1003用药时间及疗效汇总20101015.xlsx')

# meta1003 <- readxl::read_excel('D:/kintor/cn_1003_metainfor.xlsx')
meta1003 <- read_delim('D:/kintor/cn_1003_meta.txt', delim = '\t')

snv_1003 <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
                               sheet = 'DNA_SNV'
                               )
sort(table(snv_1003$GeneSymbol))

qc_1003 <- readxl::read_excel('D:/kintor/Kintor-GT0918-CN-1003_Data_03292019.xlsx',
                              sheet = 'cfDNA_QC'
                              )
cnv_1003 <- readxl::read_excel('D:/kintor/Kintor-GT0918-CN-1003_Data_03292019.xlsx',
                               sheet = 'cfDNA_CNV'
                               ) %>% 
  separate(col = SampleID, into = c('PatientID','StudyVisit'), sep = '_', remove = F)

snv_1003 %<>% separate(col = SampleID, into = c('PatientID','StudyVisit'), sep = '_', remove = F)
cfrna_splicing <- readxl::read_excel('D:/kintor/Kintor-GT0918-CN-1003_Data_03292019.xlsx',
                                     sheet = 'cfRNA_Splicing')

df1003 <- snv_1003 %>% left_join(meta1003, by = c('PatientID'='受试者编号')) %>% 
  rename('time'=`time(days)`)
df1003 %<>% mutate(time_d = ifelse(StudyVisit=='C1D1', time-0, ifelse(
  StudyVisit=='C4D1', time-84, ifelse(
    StudyVisit=='C7D1', time-84*2, ifelse(
      StudyVisit=='C10D1', time-84*3, ifelse(
        StudyVisit=='C13D1', time-84*4, ifelse(
          StudyVisit=='C16D1', time-84*5, ifelse(
            StudyVisit=='EOT', 0, NA
          )
        )
      )
    )
  )
))) %>% arrange(time_d) %>% 
  mutate(time_d=ifelse(time_d<0,0,time_d))


# 用药期长短 -------------------------------------------------------------------
# 统计测序样本数目


bar <- df %>% select(c(SubjectID, ExternalID, StudyVisit, Variant_Value, 
                       Hugo_Symbol, HGVSc, HGVSp, time, stat)) %>% 
  arrange(desc(time)) %>% 
  unite('variant', c(Hugo_Symbol, HGVSp), sep = ':', remove = T)

t <- bar %>% select(SubjectID, time, stat) %>% distinct() #按照t的分布情况划分
shapiro.test(t$time)
# 
long <- t$SubjectID[1:15]
short <- t$SubjectID[16:44]


# TP53 mutation -----------------------------------------------------------

gene_mutation_count <- function(cut, gene){
  df1003 %>% filter(PatientID %in% cut) %>% filter(GeneSymbol==gene) %>% 
    group_by(HGVS) %>% 
    # count()
    select(SampleID, Clinvar) %>% arrange(SampleID)
}
count_sample_gene <- function(cut, gene){
  df1003 %>% filter(PatientID %in% cut) %>% filter(GeneSymbol==gene) %>% 
    group_by(PatientID) %>% count()
}


t <- df1003 %>% select(PatientID, `time`) %>% arrange(`time`) %>% distinct()
# short <- c(t$PatientID[1:9],c('08001', '08003', '03001'))
short <-t$PatientID[1:6]
long <- t$PatientID[7:15]

aa <- gene_mutation_count(short, 'TP53')
count_sample_gene(short, 'TP53')

write_delim(aa, 'D:/kintor/mutation2.txt', delim = '\t')

# CNV analysis
df_2 <- cnv_1003 %>% left_join(meta1003, by = c('PatientID'='受试者编号'))
df_2 %>% filter(PatientID  %in% short) %>% filter(GeneSymbol =='AR') %>% 
  group_by(PatientID ) %>% count()

df2 %>% filter(PatientID  %in% short) %>% filter(GeneSymbol =='AR') %>% 
  group_by(`Variant_Classification`) %>% 
  # count()
  select(ExternalID) %>% arrange(ExternalID)


# AR splcing
splicing_1003 <- cnv_1003 <- readxl::read_excel('D:/kintor/Kintor-GT0918-CN-1003_Data_03292019.xlsx',
                                                sheet = 'cfRNA_Splicing'
) %>% 
  separate(col = SampleID, into = c('PatientID','StudyVisit'), sep = '_', remove = F)
df_3 <- splicing_1003 %>% left_join(meta1003, by = c('PatientID'='受试者编号'))



# survival analysis -------------------------------------------------------
fit <- survfit(Surv(`time(days)`, statues) ~ `剂量`, data = meta1003)
summary(fit)
ggsurvplot(fit, pval = T)

# TP53
bb <- c(
  'c.673-2A>G',
  'c.737T>G',
  'c.524G>A',
  'c.839G>T',
  'c.533delinsCA',
  'c.527G>T',
  'c.527G>C',
  'c.473G>A'
) # 出现bb的为预后较差的
bb <- c("TP53:c.673-2A>G",
        "TP53:p.Met246Arg",
        "TP53:p.Arg175His",
        "TP53:p.Arg280Lys",
        "TP53:p.H178Pfs*2",
        "TP53:p.Cys176Phe",
        "TP53:p.Cys176Ser",
        "TP53:p.Arg158His"
)


# 对于PatientID而言，只要患者出现yes的TP53突变，不论时期，便定义为yes
# 对于时期的考虑，从两方面出发，一是按照SampleID进行，二是只分析入组时的状态
df_p53 <- df1003 %>% mutate(stat_p53 = if_else((GeneSymbol =='TP53') & 
                                                 (HGVS %in% bb), 'yes', 'no')) %>% 
  select(SampleID, PatientID, HGVSc, `time`, stat_p53, statues, time_d) %>% 
  arrange(SampleID) %>% 
  select(PatientID, stat_p53) %>% distinct() %>% 
  group_by(PatientID) %>% 
  top_n(1)

d_tp53 <- meta1003 %>% right_join(df_p53, by = c('受试者编号'='PatientID'))

fit <- survfit(Surv(`time(days)`, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

# 按照SampleID
df_p53 <- df1003 %>% mutate(stat_p53 = if_else((GeneSymbol =='TP53') & 
                                                 (HGVS %in% bb), 'yes', 'no')) %>% 
  select(SampleID, PatientID, HGVS, `time`, stat_p53, statues, time_d) %>% 
  arrange(SampleID) %>% distinct() %>% 
  select(SampleID, PatientID, time_d, stat_p53,statues) %>% distinct() %>% 
  select(SampleID, stat_p53) %>% distinct() %>% 
  group_by(SampleID) %>% 
  top_n(1) %>% ungroup()

foo <- df1003 %>% select(SampleID, time, time_d, statues) %>% distinct()
d_tp53 <- df_p53 %>% left_join(foo, by = c('SampleID'))
d_tp53 %<>% filter(time_d > 0) # 似乎应该删去0值的存在
fit <- survfit(Surv(time_d, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)


# 只考虑患者入组时的状态
ruzu <- df1003 %>% filter(StudyVisit=='C1D1')
df_p53 <- ruzu %>% mutate(stat_p53 = if_else((GeneSymbol =='TP53') & 
                                               (HGVS %in% bb), 'yes', 'no')) %>% 
  select(SampleID, PatientID, HGVS, `time`, stat_p53, statues, time_d) %>% 
  arrange(SampleID) %>% 
  select(PatientID, stat_p53) %>% distinct() %>% 
  group_by(PatientID) %>% 
  top_n(1)
d_tp53 <- meta1003 %>% right_join(df_p53, by = c('受试者编号'='PatientID'))

fit <- survfit(Surv(`time(days)`, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)




# AR mutation

dd <- c()








