#
library(survival)
library(survminer)
library(tidyverse)
library(lubridate)
library(magrittr)
library(ComplexHeatmap)


# read file --------------------------------------------------------------

# cfdna_snv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
#                                 sheet = 'DNA_SNV')
# 
# cfdna_cnv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
#                                 sheet = 'DNA_CNV')
# cfdna_qc <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
#                                sheet = 'DNA_QC_Summary')
# cfrna_splicing <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
#                                      sheet = 'RNA_Splicing')
# cfrna_snv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_data_20210520.xlsx',
#                                 sheet = 'RNA_SNV_Confirm')

cfdna_snv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                sheet = 'DNA_SNV')

cfdna_cnv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                sheet = 'DNA_CNV')
cfdna_qc <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                               sheet = 'DNA_QC_Summary')
cfrna_splicing <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                     sheet = 'RNA_Splicing')
cfrna_snv <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                sheet = 'RNA_SNV_Confirm')
psa <- read_delim('D:/cn1003/z1.txt', delim = '\t')
drug_infor <- readxl::read_excel('D:/cn1003/1003_drug_history.xlsx',
                                 sheet = 'Sheet1') %>% 
  select(Subject, drug)

metadata <- read_delim('D:/cn1003/cn_1003_meta.txt', delim = '\t') %>% 
  left_join(drug_infor, by = c('PatientID'='Subject')) %>% 
  arrange(time)
Hmisc::describe(metadata$time)
describer::describe(metadata$time)
Rmisc::CI(metadata$time,
          ci = 0.95
)
ggplot(metadata, aes(y=time)) + geom_boxplot() + theme_bw()
ggplot(metadata, aes(x=time)) + geom_histogram() + theme_bw()

# write_delim(d, 'D:/kintor/tidy_metadata.txt', delim = '\t')

df <- cfdna_snv %>% left_join(metadata, by = c('SubjectID'='PatientID')) %>% 
  mutate(StudyVisit=if_else(StudyVisit=='End of Treatment', 'EOT', StudyVisit))

df %<>% mutate(time_d = ifelse(StudyVisit=='C1D1', time-0, ifelse(
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

# overall mutation, 统计各个mutation的个数分布
df %>% unite(col = 'mut', c(GeneSymbol,HGVSp_Short), remove = F) %>% 
  group_by(mut) %>% count() %>% arrange(desc(n)) %>% filter(n>=2)

# QC
table(cfdna_qc$LibraryQC)
table(cfdna_qc$Proceeded)
table(cfdna_qc$NGSQC)
cfdna_qc %>% filter(NGSQC %in% c('FAIL','AT RISK'))
cfdna_qc %>% filter(is.na(TotalSeqReads))


#
# suvival analysis --------------------------------------------------------

df_for_cfdna <- df %>% select(c(SubjectID, ExternalID, StudyVisit, GeneSymbol,
                                VariantFreq, HGVSc, HGVSp_Short, dosing, time, statues
)) %>% 
  distinct() %>% 
  unite('variant', c(GeneSymbol, HGVSp_Short), sep = ':', remove = F) %>% 
  # group_by(variant) %>% filter(n()>2) %>% 
  mutate(dosing = ifelse(dosing=='100', 1, 2))

df_for_all <- metadata %>% 
  mutate(dosing = ifelse(dosing=='100', 1, if_else(
    dosing=='200', 2,3
  )))

df_for_all$dosing <- factor(df_for_all$dosing, levels = c(1,2,3), 
                            labels = c('100mg', '200mg', '300mg'))

# km plot
metadata1 <- metadata %>% filter(!PatientID %in% first_id)
fit <- survfit(Surv(time, statues) ~ dosing, data = metadata1)
summary(fit)
surv_summary(fit)
ggsurvplot(fit, pval = T)
ggsurvplot(fit, 
           data = metadata1,
           surv.median.line = "hv",
           legend.title = "class",
           # legend.labs = c("Male", "Female"),
           pval = TRUE,
           pval.method = TRUE,
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           # palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_survminer() +
             theme(legend.text = element_text(size = 15),
                   legend.title = element_text(size = 15)
             )
)

# 不同批次
metadata1 <- metadata %>% mutate(batch = if_else(PatientID %in% first_id,
                                                 'first', 'second'
                                                 )) %>% 
  arrange(time)
fit <- survfit(Surv(time, statues) ~ batch, data = metadata1)
metadata1$PatientID <- factor(metadata1$PatientID, levels = metadata1$PatientID)
ggplot(metadata1, aes(x = batch, y = time, color=batch)) + geom_boxplot() +
  theme_bw() +
  ggsci::scale_color_aaas()

# different drug
metadata1 <- metadata %>% filter(!drug %in% c('NA','C+B'))
fit <- survfit(Surv(time, statues) ~ drug, data = metadata1)

ggsurvplot(fit, 
           data = metadata1,
           surv.median.line = "hv",
           legend.title = "class",
           # legend.labs = c("Male", "Female"),
           pval = TRUE,
           pval.method = TRUE,
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = 'jama',
           ggtheme = theme_survminer() +
             theme(legend.text = element_text(size = 15),
                   legend.title = element_text(size = 15)
             ))

# cox
reg <- coxph(Surv(time, statues) ~ dosing, 
             data = df_for_all)
summary(reg)
survminer::ggforest(reg, data = df_for_all)
cox.zph(reg)


# 设置用药期长??? ----------------------------------------------------------
foo <- df %>% select(c(SubjectID, ExternalID, StudyVisit, VariantFreq, 
                       GeneSymbol, HGVSc, HGVSp_Short, time, statues)) %>% 
  arrange(desc(time)) %>% 
  unite('variant', c(GeneSymbol, HGVSp_Short), sep = ':', remove = T)

t <- foo %>% select(SubjectID, time, statues) %>% distinct() #按照t的分布情况划???
shapiro.test(t$time)
Hmisc::describe(t$time)
# 
long <- t$SubjectID[1:19]
short <- t$SubjectID[20:38]



# tp53 --------------------------------------------------------------------

# 统计long/short中出现TP53的患者数
count_sample_gene <- function(cut, gene){
  df %>% filter(SubjectID %in% cut) %>% filter(GeneSymbol==gene) %>% 
  group_by(SubjectID) %>% count()
  }
count_sample_gene(long, 'TP53')

cout_ExternalID_gene <- function(df, gene, cut){
  df %>% filter(GeneSymbol==gene) %>% 
    filter(SubjectID %in% cut) %>% 
    group_by(ExternalID) %>% 
    summarise(n = n())
}
cout_ExternalID_gene(df, 'TP53', long)
# 统计出现某基因突变的全部样本
df %>% filter(GeneSymbol=='AR') %>% 
  filter(HGVSp_Short=='p.W742L') %>% 
  select(SubjectID) %>% 
  distinct()


# 统计long/short中tp53突变情况
gene_mutation_count <- function(cut, gene){
  df %>% filter(SubjectID %in% cut) %>% filter(GeneSymbol==gene) %>% 
  # group_by(HGVSp_Short) %>%
  # group_by(HGVSc) %>%
  # count()
  select(HGVSp_Short,ExternalID, Clinvar, CLIN_SIG, time, time_d, statues,
         AltDepth, VariantFreq,Mutation_Status) %>% 
    # filter(VariantFreq>1) %>%
    arrange(ExternalID) %>% 
    distinct()
  }
bar <- gene_mutation_count(short, 'TP53') %>% data.frame() %>% 
  distinct() %>% 
  arrange(time)

write_delim(bar, 'D:/cn1003/mutation.txt', delim = '\t')

# count gene
tp53_cnv_id <- df2 %>% 
                filter(Hugo_Symbol=='TP53') %>% 
                select(SubjectID) %>% distinct() %>% pull(SubjectID)

tp53_mutation_id <- df  %>% 
                     filter(GeneSymbol=='TP53') %>% 
                     filter(VariantFreq>0) %>% 
                     select(SubjectID) %>% distinct() %>% pull(SubjectID)
length(unique(c(tp53_cnv_id, tp53_mutation_id)))


# TP53 野生型和突变型 生存曲线(C1D1期)
c1d1.samples <- df %>% filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, VariantFreq, GeneSymbol,HGVSp_Short, time, statues) %>% 
  filter(GeneSymbol=='TP53') %>% select(SubjectID) %>% distinct()

m1 <- metadata %>% mutate(tp53=if_else(PatientID %in% c1d1.samples$SubjectID,
                                 'mutant','wildtype')) %>%
  select(PatientID, time, statues, tp53)
  filter(!PatientID %in% c('02003','14002'))
# m1 <- df %>% select(SubjectID, time, statues) %>% distinct() %>%
#   mutate(tp53=if_else(SubjectID %in% c1d1.samples$SubjectID, 'mutant','wild'))
fit <- survfit(Surv(time, statues) ~ tp53, data = m1)
ggsurvplot(fit, pval = T)

# TP53 野生型和突变型 生存曲线(C1D1期), 只考虑pathogenic, likly pathognic
c1d1.samples <- df %>% 
  # filter(StudyVisit=='C1D1') %>%
  select(SubjectID, VariantFreq, GeneSymbol,HGVSp_Short, Clinvar,CLIN_SIG,
         time, statues) %>% 
  filter(GeneSymbol=='TP53') %>% 
  filter(Clinvar %in% c('Pathogenic','Likely pathogenic',
                        'Pathogenic/Likely_pathogenic')) %>% 
  select(SubjectID) %>% distinct()



#
# AR ----------------------------------------------------------------------

count_sample_gene(long, 'AR')

# 统计long/short中AR突变情况
cout_ExternalID_gene(df, 'AR', long)

# 统计AR 各mutation的突变情况, 在多少sample中检测到以及在C1期检测到
df %>% filter(GeneSymbol=='AR') %>% 
  select(SubjectID, ExternalID, StudyVisit,VariantFreq, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short, StudyVisit) %>% 
  summarise(n = n()) %>% arrange(n) #共24个sample 
df %>% filter(GeneSymbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, ExternalID, VariantFreq, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short) %>% 
  summarise(n = n()) %>% arrange(n) #15个subjectID
# 统计所有突变 在 不同采样期出现的数目
aaaa <- df %>% filter(GeneSymbol=='AR') %>% 
  select(SubjectID, ExternalID, StudyVisit,VariantFreq, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short, StudyVisit) %>% 
  summarise(n = n()) %>% arrange(n) %>% 
  pivot_wider(names_from = StudyVisit, values_from = n)
write_delim(aaaa, 'D:/cn1003/genemutation.txt', delim = '\t')

# AR 野生型和突变型 生存曲线(C1D1期)
# 考虑 W742C为 获益 突变
c1d1.samples <- df %>% filter(StudyVisit=='C1D1') %>%
  filter(GeneSymbol=='AR') %>% 
  # filter(HGVSp_Short != 'p.W742C') %>%
  select(SubjectID, time) %>% distinct() %>% 
  pull(SubjectID)
m1 <- metadata %>% mutate(AR=if_else(PatientID %in% c(c1d1.samples), 
                                       'mutant','wildtype')) %>% 
  select(PatientID, time, statues, AR)
  filter(PatientID %in% first_id) %>% 
  filter(!PatientID %in% c('02003','14002'))
table(m1$AR)
m1 <- m1 %>% mutate(AR=ifelse(PatientID %in% c('02003','14002'), 'wildtype', AR))

fit <- survfit(Surv(time, statues) ~ AR, data = m1)
ggsurvplot(fit, pval = T)


# 对于ARcnv而言-----------------------------------------------
sort(table(cfdna_cnv$Hugo_Symbol)) #发生CNV最多次的为AR,PIK3CA,PB1,PTEN,AKT3,MYC
df2 <- cfdna_cnv %>% left_join(metadata, by = c('SubjectID'='PatientID'))
# long/short中CNV出现数目
df2 %>% filter(SubjectID %in% short) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(SubjectID) %>% count()

# sample样本，每个样本发生CNV的个数
df2 %>% group_by(ExternalID) %>% 
  count() %>% 
  arrange(desc(n))
# 某基因发生CNV的所有样本
hxd2 <- df2 %>% filter(Hugo_Symbol=='AR')


# 分时期对某些基因CNV 进行展示
i_gene <- Hmisc::Cs(AR,PIK3CA,RB1,PTEN,AKT3,MYC,BRCA2,POLE,BRAF)
i_cnv <- df2 %>% filter(StudyVisit=='C1D1') %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(time) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))
# i_cnv <- i_cnv[, order(names(i_cnv))]
ff <- metadata %>% select(PatientID,time) %>% 
  arrange(time) %>% 
  left_join(i_cnv, by = c('PatientID'='SubjectID')) %>% 
  filter(!PatientID %in% first_id)
M <- ff %>% 
  select(-time) %>% 
  column_to_rownames('PatientID')
ha <- columnAnnotation(bar = anno_barplot(ff$time))


# 所有发生AR 的样本
aa <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  select(ExternalID, Variant_Value, time, statues) %>% 
  arrange(time) %>% filter(ExternalID != '14009_EOT')
write_delim(aa, 'D:/cn1003/mutation1.txt', delim = '\t')
# M <- aa %>% select(-c(time, statues)) %>% 
#   column_to_rownames('ExternalID')
# ha <- rowAnnotation(foo = anno_barplot(aa$time))


# CNV survival analysis
cnv_gene<- c('RB1', 'AR', 'BRAF','BRCA2')

cnv.s <- function(gene){
  c1d1 <- df2 %>% filter(StudyVisit=='C1D1') %>% 
    filter(Hugo_Symbol %in% gene) %>% 
    select(SubjectID) %>% distinct()
  c1d1$SubjectID
}

cnv_id <- unique(cnv.s('AR'))
m1 <- metadata %>% select(PatientID,time, statues) %>% 
  distinct() %>% 
  # filter(PatientID %in% first_id) %>% 
  mutate(gene=ifelse((PatientID %in% cnv_id),'mutant', 'wildtype'))
  
fit <- survfit(Surv(time, statues) ~ gene, data = m1)
ggsurvplot(fit, pval = T)


# 选择的CNV包括AR,RB1,PIK3CA,PTEN,AKT3,BRAF,BRCA2
cnv.RB1 <- cnv.s('RB1')
cnv.PIK3CA <- cnv.s('PIK3CA')
cnv.AR <- cnv.s('AR')
cnv.BRCA2 <- cnv.s('BRCA2')
cnv.PTEN <- cnv.s('PTEN')
cnv.AKT3 <- cnv.s('AKT3')
cnv.PIK3CB <- cnv.s('PIK3CB')

m1 <- metadata %>% select(PatientID, time, statues) %>% distinct() %>%
  mutate(RB1=if_else(PatientID %in% cnv.RB1, 1,0)) %>%
  mutate(PIK3CA=if_else(PatientID %in% cnv.PIK3CA, 1,0)) %>% 
  mutate(BRCA2=if_else(PatientID %in% cnv.BRCA2, 1,0)) %>% 
  mutate(AR=if_else(PatientID %in% cnv.AR, 1,0)) %>% 
  mutate(PTEN=if_else(PatientID %in% cnv.PTEN, 1,0)) %>% 
  mutate(AKT3=if_else(PatientID %in% cnv.AKT3, 1,0)) %>% 
  mutate(PIK3CB=if_else(PatientID %in% cnv.PIK3CB, 1,0))


reg <- coxph(Surv(time, statues) ~ AR,
             data = m1)
summary(reg)
survminer::ggforest(reg, data = m1,
                    fontsize = 1.2
                    )


#
# AR Splicing-----------------------------------
# 对于AR splicing而言，只对发生splicing的样本进行分
df3 <- cfrna_splicing %>% left_join(metadata, by = c('SubjectID'='PatientID'))
sort(table(df3$NameAST)) #AR-2_3,AR-3_4
# long/short中出现splicing的频
df3 %>% filter(SubjectID %in% long) %>% 
  group_by(SubjectID) %>% count()
# 某基因出现的样本数目
hxd3 <- df3 %>% filter(NameAST == 'AR-V7')


# AR splicing的分析以AR-V7为主
# 首先还是出现的频???
df3 %>% filter(NameAST == 'AR-V7') %>% 
  filter(SubjectID %in% long) %>% 
  distinct(SubjectID)
# 发生ARV7的样本，其的变异情况
arv7 <- df3 %>% 
  filter(NameAST %in% c('AR-V7'
                        # 'AR-V3', 'AR-V9'
                        )) %>% 
  filter(SupportReads >= 3) %>% 
  select(SubjectID, ExternalID, StudyVisit, 
         NameAST, SupportReads, time, statues) %>% 
  arrange(desc(time))

write_delim(arv7, 'D:/cn1003/arv7.txt', delim = '\t')
ar_filter <- arv7 %>% filter(time > 300)
aaid <- ar_filter$ExternalID

(df_arv7 <- df %>% filter(ExternalID %in% aaid) %>% 
  select(SubjectID, ExternalID, StudyVisit, GeneSymbol, time) %>% 
  arrange(SubjectID))
df_arv7 %>% data.frame() %>% filter(GeneSymbol=='TP53')


# plot  -------------------------------------------------------------------

# barplot
aa <- metadata %>% select(PatientID, time, drug) %>% distinct() %>% 
  mutate(still=if_else(PatientID %in% first_id, 'first', 'second')) %>% 
  arrange(time)
still_drug <- c('14007','05008','14009','02003','14002','01009')

aa$PatientID <- factor(aa$PatientID, levels = aa$PatientID)
# aa %<>% mutate(still=if_else(PatientID %in% first_id, 'first', 'second'))

ggplot(aa) + geom_bar(aes(x = PatientID, y = time), 
                      stat = 'identity')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = c(0.1,0.9),
        legend.direction = 'horizontal'
        ) + ggsci::scale_fill_aaas() +
  geom_hline(yintercept = 372, linetype = 2)

# heatmap
library(ComplexHeatmap)

# 挑选样本某个基因最高突变频率的突变
df %>% unite('mut',c(GeneSymbol, HGVSp_Short), remove = F) %>% 
  filter(GeneSymbol=='TP53') %>% 
  select(SubjectID, VariantFreq, mut) %>% 
  group_by(SubjectID) %>% 
  slice_max(VariantFreq) %>% 
  pivot_wider(names_from = mut, values_from = VariantFreq)
  
# C1D1期某个基因的突变频率
filter_gene <- function(gene, visit='C1D1'){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, VariantFreq) %>% 
    select(SubjectID,VariantFreq) %>% 
    group_by(SubjectID) %>% 
    slice_max(VariantFreq) %>% 
    dplyr::rename(!!sym(gene) := VariantFreq)
  res
  }
mutation <- filter_gene('TP53')
mutation_ar <- filter_gene('AR')
mutation_atm <- filter_gene('ATM')
mutation_pik3ca <- filter_gene('PIK3CA')

ar_cnv <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')

metadata %<>% left_join(psa, by = c('PatientID'='SUBJID')) %>% 
  arrange(max) %>% 
  mutate(max = if_else(max>1,1,max))
final_data <- metadata %>% select(PatientID,max) %>% 
  left_join(mutation, by = c('PatientID'='SubjectID')) %>% 
  left_join(mutation_atm, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_pik3ca, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_ar, by = c('PatientID'='SubjectID')) %>%
  left_join(ar_cnv, by = c('PatientID'='SubjectID'))

M <- final_data %>% 
  select(-max) %>%
  column_to_rownames('PatientID')

ha <- columnAnnotation(bar = anno_barplot(metadata$max))

# 下面为某基因全部突变
filter_aa <- function(gene, aa){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    filter(StudyVisit=='C1D1') %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, VariantFreq) %>% 
    select(SubjectID,VariantFreq) %>% 
    group_by(SubjectID) %>% 
    slice_max(VariantFreq) %>% 
    rename((!!aa) := VariantFreq)
  res
}
ar_w742c <- filter_aa('AR','p.W742C')

ar_gene <- c('p.W742C', 'p.W742L', 'p.H875Y','p.T878A','p.L702H'
             # 'p.I883M', 'p.P893S','p.T878S','p.V716M','p.S889G'
             )
final_data <- metadata %>% select(PatientID,max)
for(i in ar_gene){
  foo <- filter_aa('AR', i)
  final_data <- final_data %>%
    left_join(foo, by = c('PatientID'='SubjectID'))
  
}


Heatmap(t(M),
        name = 'ratio',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        # col = circlize::colorRamp2(c(0, 2, 4), c("navy", "white", "firebrick3")),
        col = circlize::colorRamp2(c(0,10), c("white", "firebrick3")),
        width = ncol(M)*unit(40, "mm"),
        height = nrow(M)*unit(1, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = F,
        # column_names_gp = gpar(fontsize=1)
)

pheatmap::pheatmap(M)

# 按照ExternalID 进行
filter_gene_e <- function(gene){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    arrange(time) %>% 
    select(ExternalID,VariantFreq) %>% 
    group_by(ExternalID) %>% 
    slice_max(VariantFreq) %>% 
    dplyr::rename(!!sym(gene) := VariantFreq)
  res
}
tp53.e <- filter_gene_e('TP53')
f1 <- df %>% select(ExternalID,time) %>% 
  distinct() %>% 
  arrange(time)
f <- f1 %>% left_join(tp53.e, by = c('ExternalID'))
M <- f %>% 
  select(-time) %>%
  column_to_rownames('ExternalID')
ha <- columnAnnotation(bar = anno_barplot(f1$time))




# 只从突变频率展开分析 --------------------------------------------------------------
# 假设高于5%的somatic突变频率为有害因素
# 只考虑C1D1
df %>% filter(VariantFreq>=5) %>% 
  filter(StudyVisit=='C1D1') %>% 
  filter(Mutation_Status=='Likely Somatic') %>% 
  arrange(SubjectID, VariantFreq) %>% 
  select(SubjectID,VariantFreq,time, statues) %>% 
  group_by(SubjectID) %>% 
  slice_max(VariantFreq)



# survival by mutation ----------------------------------------------------

aa <- as.character(quote(
  c(p.V143M,
    p.R280I,
    p.C176F,
    p.C176S,
    p.K139N,
    p.F270V,
    p.S215I,
    p.M246R,
    p.R158H,
    p.Y220D,
    p.N288D,
    p.D281N
  )
))[-1]
aa <- c(aa, 'p.H178Pfs*3','p.M133Ifs*37')

# 按照SubjectID
df_p53 <- df %>% mutate(stat_p53 = ifelse(GeneSymbol=='TP53' & 
                                            HGVSp_Short %in% aa,
                                          'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)


d_tp53 <- metadata %>% right_join(df_p53, by = c('PatientID'='SubjectID'))

fit <- survfit(Surv(time, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

# 按照ExternalID
foo <- df %>% select(ExternalID, time, time_d, statues) %>% distinct()
df_p53 <- df %>% mutate(stat_p53 = ifelse(GeneSymbol=='TP53' & 
                                            HGVSp_Short %in% aa,
                                          'yes', 'no')) %>%
  # select(SubjectID, stat_p53, ExternalID, time, time_d, stat) %>%
  select(ExternalID, stat_p53) %>%
  distinct() %>% 
  arrange(ExternalID, desc(stat_p53)) %>% 
  group_by(ExternalID) %>% 
  top_n(1) %>% ungroup()
d_tp53 <- df_p53 %>% left_join(foo, by = c('ExternalID'))
fit <- survfit(Surv(time_d, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

write_delim(d_tp53, 'D:/cn1003/d_tp53.txt', delim = '\t')
# d_tp53 <- read_delim('D:/cn1003/d_tp53.txt', delim = '\t')


# 只考虑患者入组时的状态
ruzu <- df %>% filter(StudyVisit=='C1D1')
df_p53 <- ruzu %>% mutate(stat_p53 = ifelse(GeneSymbol=='TP53' & 
                                              HGVSp_Short %in% aa,
                                            'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)

d_tp53 <- metadata %>% right_join(df_p53, by = c('PatientID'='SubjectID'))
fit <- survfit(Surv(time, statues) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)
## 考虑突变频率，设没有TP53突变的突变频率为0.
ruzu %>% filter(GeneSymbol=='TP53') %>% select(SubjectID, VariantFreq, HGVSp_Short)


# 按照入组时候的状态，TP53、AR mutation、AR CNV、共同进行cox回归
ruzu <- df %>% filter(StudyVisit=='C1D1')
ruzu2 <- df2 %>% filter(StudyVisit=='C1D1') %>% 
  filter(Hugo_Symbol == 'AR')
ruzu3 <- df3 %>% filter(StudyVisit=='C1D1') %>% 
  filter(NameAST %in% c('AR-V7','AR-V3','AR-V9'))
df_all <- ruzu %>% mutate(stat_p53 = ifelse(GeneSymbol=='TP53' & 
                                              HGVSp_Short %in% aa,
                                            'yes', 'no'),
                          stat_AR = if_else(GeneSymbol=='AR' &
                                              HGVSp_Short %in% c('p.W742C'),
                                            'yes','no'
                                            )
                          ) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)


#
# re-analysis of cfdna_SNV ------------------------------------------------

cal_gene_appear <- function(gene){
  foo <- df %>% filter(GeneSymbol==gene) %>% group_by(ExternalID) %>% 
    count()
  return(dim(foo)[1])
}
cal_gene_appear('TP53')

y <- vector()
for (i in unique(df$GeneSymbol)) {
  aa <- cal_gene_appear(i)
  names(aa) <- i
  print(aa)
  y <- c(y, aa)
}

# 统计发生突变的ExternalID
count_ExternalID <- function(df, gene){
  foo <- df %>% filter(GeneSymbol==gene) %>% group_by(ExternalID) %>% 
    summarise(n = n())
  return(dim(foo)[1])
}

count_ExternalID(cfrna_snv, 'TP53')

x <- vector()
for (i in unique(cfrna_snv$GeneSymbol)) {
  aa <- count_ExternalID(cfrna_snv, i)
  names(aa) <- i
  print(aa)
  x <- c(x, aa)
}
sort(x)



# RNA analysis ------------------------------------------------------------

m1 <- metadata %>% select(PatientID, time)
cfrna_exp <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                sheet = 'RNA_Expression') %>% 
  left_join(m1, by = c('SubjectID'='PatientID')) %>% 
  arrange(time)

# 对ExternalID 即samples名称进行
fpkm <- cfrna_exp %>% select(-c(SubjectID,StudyVisit,RequisitionID, time)) %>% 
  column_to_rownames('ExternalID') %>% t() %>% 
  as.data.frame()
# fpkm <- na.omit(fpkm)
fpkm %<>% replace(is.na(.), 0)
cc <- cfrna_exp %>% mutate(class = if_else(SubjectID %in% short, 'short', 'long'))

anno <- cc %>% select(ExternalID, class) %>% 
  column_to_rownames('ExternalID')

p_heat <- pheatmap::pheatmap(
  log10(fpkm + 1),
  border = F,
  cluster_cols = F,
  cluster_rows = T,
  show_rownames = T,
  show_colnames = T,
  fontsize_col = 6,
  # annotation_col = anno
)
# use complexheatmap
col_anno <- cfrna_exp %>% select(ExternalID, time) %>% 
  column_to_rownames('ExternalID')
ha <- columnAnnotation(time = anno_barplot(col_anno))

Heatmap(
  log10(fpkm + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=1),
  cluster_columns = F,
  cluster_rows = T
)


# only for C1D1
fpkm_c1d1 <- cfrna_exp %>% filter(StudyVisit=='C1D1') %>% 
  select(-c(ExternalID,StudyVisit,RequisitionID,time)) %>% 
  column_to_rownames('SubjectID') %>% t() %>% 
  as.data.frame()
fpkm_c1d1 %<>% replace(is.na(.), 0)

# anno <- cc %>% filter(StudyVisit=='C1D1') %>%
#   select(SubjectID, class) %>% 
#   column_to_rownames('SubjectID')
# p_heat <- pheatmap::pheatmap(
#   log10(fpkm_c1d1 + 1),
#   cluster_cols = T,
#   cluster_rows = F,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = anno
# )
# use complexheatmap
col_anno <- cfrna_exp %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, time) %>% 
  distinct() %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(time = anno_barplot(col_anno))
Heatmap(
  log10(fpkm_c1d1 + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=1),
  cluster_rows = T,
  cluster_columns = F
)


# 只选择感兴趣的基因
i_gene <- Hmisc::Cs(TP53,AR,PIK3CA,ATM,SPOP)
s <- cfrna_exp %>% 
  arrange(time) %>% 
  filter(StudyVisit=='C1D1') %>%
  select(any_of(i_gene), SubjectID) %>% 
  column_to_rownames('SubjectID') %>% t() %>% 
  as.data.frame()
col_anno <- cc %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, time) %>% 
  distinct() %>% 
  arrange(time) %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(time = anno_barplot(col_anno))
Heatmap(
  log10(s + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=10),
  cluster_rows = T,
  cluster_columns = F
)


# survival, 以中位值做区分
df4 <- cfrna_exp %>% left_join(d, by = c('SubjectID'='pt. ID')) %>% 
  mutate(status = if_else(stat==2,1,0)) %>% 
  arrange(TP53)

reg <- coxph(Surv(time, stat) ~ TP53, 
             data = df4)
summary(reg)
survminer::ggforest(reg, data = df4)
# 寻找TP53 best cutoff
s.cut <- surv_cutpoint(
  data = df4,
  time = 'time',
  event = 'status',
  variables = 'TP53'
)
summary(s.cut)
plot(s.cut, "TP53", palette = "npg")



# cfrna SNV-------------------------------------------------------------------

cfrna_snv_others <- readxl::read_excel('D:/cn1003/GT0918_CN_1003_cfDNA_cfRNA_cumulative_data_20210521.xlsx',
                                sheet = 'RNA_SNV_Others')

# cfrna confirm 最多的为TP53
cfrna_snv %>% group_by(GeneSymbol) %>% 
  count() %>% 
  arrange(desc(n))

correlation::cor_test(cfrna_snv, 'VarFreqDNA', 'VarFreqRNA',
                      method = 'pearson'
                      ) # 一致的突变频率相关性一般，且RNA的有许多为0，

ggplot(cfrna_snv,  aes(x = VarFreqDNA, y = VarFreqRNA)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

ggstatsplot::ggscatterstats(cfrna_snv, x = VarFreqDNA, y = VarFreqRNA)

# 进一步验证三文件突变位点的一致性
rna_key <- (cfrna_snv %>% 
  unite('key', c(ExternalID,GeneSymbol,HGVSp), remove = F))$key %>% 
  unique()

dna_key <- (df %>% 
  unite('key', c(ExternalID,GeneSymbol,HGVSp), remove = F))$key %>% 
  unique()

rna_other_key <- (cfrna_snv_others %>% 
  unite('key', c(ExternalID,GeneSymbol,HGVSp), remove = F))$key %>% 
  unique()

sum(rna_other_key %in% dna_key)
sum(rna_key %in% dna_key)
print(length(rna_key))
print(length(dna_key))
sum(rna_key %in% rna_other_key)


#
# 其他mutation和CNV ----------------------------------------------------------
# 全部mutation再C1D1期的分布
c1d1.all.mutation <- df %>% filter(StudyVisit == 'EOT') %>% 
  filter(VariantFreq >= 0) %>% 
  group_by(GeneSymbol) %>% 
  count() %>% arrange(desc(n))

# write_delim(c1d1.all.mutation,'D:/cn1003/c1d1.all.mutation.txt', delim = '\t')

df %>% filter(StudyVisit == 'EOT') %>% 
  filter(VariantFreq >= 0.5) %>% 
  select(SubjectID, GeneSymbol) %>% distinct() %>% 
  group_by(GeneSymbol) %>% count() %>% arrange(desc(n)) %>% 
  write_delim(.,'D:/cn1003/c1d1.all.mutation.txt', delim = '\t')


df2 %>% 
  filter(StudyVisit=='EOT') %>% 
  group_by(Hugo_Symbol) %>% 
  count() %>% arrange(desc(n)) %>% 
  write_delim(.,'D:/cn1003/c1d1.all.cnv.txt', delim = '\t')


# 某基因在C1D1期出现与否的生存分析
# AR,TP53,ATM,DNMT3A,ARID1A,BRCA2,SPOP,TET2
repair_gene <- Hmisc::Cs(BRCA2,MSH2,MSH6)

c1d1.samples <- df %>% 
  filter(StudyVisit=='C1D1') %>%
  filter(VariantFreq >= 0) %>% 
  select(SubjectID, VariantFreq, GeneSymbol,HGVSp_Short, Clinvar,
         CLIN_SIG, time, statues) %>% 
  filter(GeneSymbol %in% repair_gene) %>% 
  # filter(Clinvar %in% c('Pathogenic','Likely pathogenic',
                        # 'Pathogenic/Likely_pathogenic')) %>%
  select(SubjectID) %>% distinct() %>% 
  pull(SubjectID)
print(length(c1d1.samples))

m1 <- metadata %>% select(PatientID, time, statues) %>% distinct() %>%
  # filter(!PatientID %in% c('02003','14002')) %>%
  mutate(Repair=if_else(PatientID %in% c1d1.samples, 'mutant','wildtype'))

fit <- survfit(Surv(time, statues) ~ Repair, data = m1)
ggsurvplot(fit, pval = T)

# cox regression
c1.gene.mutate <- function(gene){
  df %>% 
    filter(StudyVisit=='C1D1') %>%
    select(SubjectID, VariantFreq, GeneSymbol,HGVSp_Short, Clinvar,
           CLIN_SIG, time, statues) %>% 
    filter(GeneSymbol==gene) %>% 
    # filter(Clinvar %in% c('Pathogenic','Likely pathogenic',
    # 'Pathogenic/Likely_pathogenic')) %>%
    select(SubjectID) %>% distinct() %>% 
    pull(SubjectID)
}

g.BRCA2 <- c1.gene.mutate('BRCA2')
g.AR <- c1.gene.mutate('AR')
g.TP53 <- c1.gene.mutate('TP53')
g.ATM <- c1.gene.mutate('ATM')
g.MSH2 <- c1.gene.mutate('MSH2')
g.MSH6 <- c1.gene.mutate('MSH6')
g.SPOP <- c1.gene.mutate('SPOP')

m1 <- metadata %>% select(PatientID, time, statues) %>% distinct() %>%
  mutate(BRCA2=if_else(PatientID %in% g.BRCA2, 1,0)) %>%
  mutate(AR=if_else(PatientID %in% g.AR, 1,0)) %>% 
  mutate(TP53=if_else(PatientID %in% g.TP53, 1,0)) %>% 
  mutate(ATM=if_else(PatientID %in% g.ATM, 1,0)) %>% 
  mutate(MSH2=if_else(PatientID %in% g.MSH2, 1,0)) %>% 
  mutate(MSH6=if_else(PatientID %in% g.MSH6, 1,0)) %>% 
  mutate(SPOP=if_else(PatientID %in% g.SPOP, 1,0))

+ATM+MSH2+MSH6
reg <- coxph(Surv(time, statues) ~ BRCA2,
             data = m1)
summary(reg)
survminer::ggforest(reg, data = m1,
                    fontsize = 1.5
                    )



# cfdna concentration -----------------------------------------------------
# 浓度变化的比???, 或者看一下浓度和time之间的相关???
# sperman correlation

haha <- cfdna_qc %>% filter(StudyVisit=='C1D1') %>% 
  arrange(SubjectID)

hehe <- t %>% arrange(SubjectID)

heha <- bind_cols(haha$DNAYield, hehe$time)
names(heha) <- c('DNA','time')
correlation::cor_test(heha, x = 'DNA', y = 'time', method = 'spearman')

# 浓度和生存的关系


# 关于mutation 突变频率



# PSA数据整理 -----------------------------------------------------------------

psa_f <- readxl::read_excel('D:/cn1003/Edata-PSA.xlsx', sheet = 'Sheet1')
m2 <- metadata %>% select(PatientID,dosing,time,statues)
m2$PatientID <- as.numeric(m2$PatientID)

patient_id <- as.numeric(metadata$PatientID)

d <- psa_f %>% dplyr::select(SUBJID, VISITNUM, LBORRES) %>% 
  filter(SUBJID %in% patient_id) %>% 
  filter(VISITNUM != 'UNPLANNED') %>% 
  # group_by(VISITNUM) %>%
  # mutate(row = row_number()) %>%
  pivot_wider(names_from = VISITNUM, values_from = LBORRES, 
              # values_fn = list(LBORRES=length)
              ) %>% 
  left_join(m2, by = c('SUBJID'='PatientID')) %>% 
  arrange(time)

d %<>% mutate(C2S=(C2D1-SCREEN)/SCREEN) %>% 
  mutate(C3S=(C3D1-SCREEN)/SCREEN) %>% 
  mutate(C4S=(C4D1-SCREEN)/SCREEN) %>% 
  mutate(C5S=(C5D1-SCREEN)/SCREEN) %>% 
  mutate(C6S=(C6D1-SCREEN)/SCREEN) %>% 
  mutate(EOTS=(EOT-SCREEN)/SCREEN) %>% 
  mutate(SUBJID=as.character(SUBJID)) %>% 
  filter(SUBJID %in% as.character(as.numeric(id_c)))
  
# write_delim(d, 'D:/cn1003/psa.txt', delim = '\t')

dd <- d %>% select(SUBJID,C2S,C3S,C4S,C5S,C6S,EOTS)

M <- dd %>% 
  column_to_rownames('SUBJID')
# M <- d %>% select(SUBJID,SCREEN,C2D1,C3D1,C4D1,C5D1,C6D1,EOT) %>% 
#   column_to_rownames('SUBJID')
ha <- columnAnnotation(time = anno_barplot(d$time))
Heatmap(t(M),
        name = 'ratio',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        col = circlize::colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
        # width = ncol(M)*unit(30, "mm"),
        height = nrow(M)*unit(5, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)

ddd <- d %>% select(C2S,C3S,C4S,C5S,C6S,EOTS) %>% 
  rowwise() %>% 
  mutate(value=max(C2S,C3S,C4S,C5S,C6S,EOTS, na.rm = T))
ddd$max <- apply(ddd, 1, function(x){x[which.max(abs(x))]})

dddd <- bind_cols(dd, (ddd %>% select(max)))
dddd$SUBJID <- factor(dddd$SUBJID, levels = dddd$SUBJID)
dddd %<>% mutate(col=if_else(max>0,'up','down'))

library(ggbreak)
library(patchwork)
ggplot(dddd, aes(x=SUBJID, y=max, fill=col)) +
  # geom_bar(stat = 'identity') + 
  geom_col() +
  scale_y_break(c(1,34), scales = 0.2, ticklabels = NULL) +
  theme_bw() +
  ggsci::scale_fill_jama() +
  ylab('PSA Ratio(%)') +
  # geom_hline(yintercept = -0.5, colour="red", linetype="dashed")
  theme(axis.text.x = element_text(angle = 90,
                                   size = 10),
        axis.text.y = element_text(size = 10)
  ) 

m2 %<>% mutate(PatientID=as.character(PatientID))
z1 <- dddd %>% left_join(m2, by = c('SUBJID'='PatientID')) %>% 
  arrange(desc(max))
z1$SUBJID <- factor(z1$SUBJID, levels = z1$SUBJID)
correlation::cor_test(z1,'max','time',method = 'spearman')

z1 %<>% mutate(t = case_when(
  time >= 365 ~ '>365',
  time < 365 ~ '<365',
  TRUE ~ 'other'
))

cnv.id <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID) %>% distinct()

z1 %<>% mutate(arcnv=if_else(
  SUBJID %in% as.numeric(cnv.id$SubjectID),'yes', 'no'
))

ggplot(z1) + geom_bar(aes(x = SUBJID, y = max, fill=arcnv), 
                      stat = 'identity')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = 'none'
  ) + ggsci::scale_fill_aaas() +
  scale_y_break(c(1,10), scales = 0.2, ticklabels = NULL) +
  ylab('Change from Baseline')

# write_delim(z1, 'D:/cn1003/z1.txt', delim = '\t')


# 重新定义EOT -----------------------------------------------------------------
# 测序样本只有38例
eot.na <- cfdna_qc %>% select(SubjectID, StudyVisit, DNAYield) %>% 
  distinct() %>% 
  pivot_wider(names_from = StudyVisit, values_from = DNAYield) %>% 
  filter(is.na(EOT))
eot.na$SubjectID

eot <- function(gene){df %>% filter(GeneSymbol==gene) %>% 
  filter(StudyVisit=='EOT') %>% 
  arrange(SubjectID, VariantFreq) %>% 
  select(SubjectID,VariantFreq) %>% 
  group_by(SubjectID) %>% 
  slice_max(VariantFreq) %>% 
  dplyr::rename(!!sym(gene) := VariantFreq)
    }

eott <- function(gene){df %>% filter(GeneSymbol==gene) %>% 
  filter(StudyVisit=='C4D1') %>% 
  arrange(SubjectID, VariantFreq) %>% 
  select(SubjectID,VariantFreq) %>%
  filter(SubjectID %in% eot.na$SubjectID) %>% 
  group_by(SubjectID) %>% 
  slice_max(VariantFreq) %>% 
    dplyr::rename(!!sym(gene) := VariantFreq)
  }

mutation <- bind_rows(eot('TP53'), eott('TP53'))
mutation_ar <- bind_rows(eot('AR'), eott('AR'))
mutation_atm <- bind_rows(eot('ATM'), eott('ATM'))
mutation_pik3ca <- bind_rows(eot('PIK3CA'), eott('PIK3CA'))

ar_cnv1 <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='EOT') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')
ar_cnv2 <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C4D1') %>% 
  select(SubjectID, Variant_Value) %>% 
  filter(SubjectID %in% eot.na$SubjectID) %>% 
  rename('AR_CNV'='Variant_Value')
ar_cnv <- bind_rows(ar_cnv1, ar_cnv2)

metadata %<>% arrange(time)
final_data <- metadata %>% select(PatientID,time) %>% 
  left_join(mutation, by = c('PatientID'='SubjectID')) %>% 
  left_join(mutation_atm, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_pik3ca, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_ar, by = c('PatientID'='SubjectID')) %>%
  left_join(ar_cnv, by = c('PatientID'='SubjectID')) %>% 
  filter(PatientID %in% id_bc)

M <- final_data %>% 
  select(-time) %>%
  column_to_rownames('PatientID')

ha <- columnAnnotation(bar = anno_barplot(final_data$time))

# PSA MUTATION
metadata %<>% left_join(psa, by = c('PatientID'='SUBJID')) %>% 
  arrange(max) %>% 
  mutate(max = if_else(max>1,1,max))
final_data <- metadata %>% select(PatientID,max) %>% 
  left_join(mutation, by = c('PatientID'='SubjectID')) %>% 
  left_join(mutation_atm, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_pik3ca, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_ar, by = c('PatientID'='SubjectID')) %>%
  left_join(ar_cnv, by = c('PatientID'='SubjectID'))

M <- final_data %>% 
  select(-max) %>%
  column_to_rownames('PatientID')

ha <- columnAnnotation(bar = anno_barplot(metadata$max))


df %>% filter(GeneSymbol=='TP53') %>% 
  arrange(SubjectID, VariantFreq) %>% 
  select(ExternalID, SubjectID, StudyVisit, VariantFreq) %>% 
  group_by(ExternalID) %>% 
  slice_max(VariantFreq) %>% 
  arrange(StudyVisit) %>% 
  pivot_wider(names_from = StudyVisit, values_from = VariantFreq) %>% 
  arrange(ExternalID)

# AR mutation
  
filter_aa <- function(gene, aa){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    filter(StudyVisit=='EOT') %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, VariantFreq) %>% 
    select(SubjectID,VariantFreq) %>% 
    group_by(SubjectID) %>% 
    slice_max(VariantFreq) %>% 
    rename((!!aa) := VariantFreq)
  res
}
filter_c4d1 <- function(gene, aa){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    filter(StudyVisit=='C4D1') %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, VariantFreq) %>% 
    select(SubjectID,VariantFreq) %>% 
    group_by(SubjectID) %>% 
    slice_max(VariantFreq) %>% 
    rename((!!aa) := VariantFreq)
  res
}

ar_gene <- c('p.W742C', 'p.W742L', 'p.H875Y','p.T878A','p.L702H')

final_data <- metadata %>% select(PatientID,time)
for(i in ar_gene){
  foo <- filter_aa('AR', i)
  final_data <- final_data %>%
    left_join(foo, by = c('PatientID'='SubjectID'))
  
}

final_data <- final_data %>% 
  filter(!PatientID %in% eot.na$SubjectID)

final_C4_EOT <- metadata %>% select(PatientID,time)
for(i in ar_gene){
  foo <- filter_c4d1('AR', i)
  final_C4_EOT <- final_C4_EOT %>%
    left_join(foo, by = c('PatientID'='SubjectID'))
  
}
final_C4_EOT <- final_C4_EOT %>% 
  filter(PatientID %in% eot.na$SubjectID)

final_data <- bind_rows(final_data, final_C4_EOT) %>% 
  arrange(time)

M <- final_data %>% 
  select(-time) %>%
  column_to_rownames('PatientID')

ha <- columnAnnotation(bar = anno_barplot(metadata$time))


# 全部的CNV 
i_gene <- Hmisc::Cs(AR,PIK3CA,RB1,PTEN,AKT3,MYC,BRCA2,POLE,BRAF)
i_cnv <- df2 %>% filter(StudyVisit=='EOT') %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(time) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))

i_cnv_2 <- df2 %>% filter(StudyVisit=='C4D1') %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(time) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  filter(SubjectID %in% eot.na$SubjectID) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))

i_cnv %<>% bind_rows(i_cnv_2)

ff <- metadata %>% select(PatientID,time) %>% 
  arrange(time) %>% 
  left_join(i_cnv, by = c('PatientID'='SubjectID')) %>% 
  filter(!PatientID %in% first_id)
M <- ff %>% 
  select(-time) %>% 
  column_to_rownames('PatientID')
ha <- columnAnnotation(bar = anno_barplot(ff$time))

Heatmap(t(M),
        name = 'ratio',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        col = circlize::colorRamp2(c(0, 2, 4), c("navy", "white", "firebrick3")),
        # col = circlize::colorRamp2(c(0,20), c("white", "firebrick3")),
        width = ncol(M)*unit(30, "mm"),
        height = nrow(M)*unit(3, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)



#
# 线性回归 --------------------------------------------------------------------

library(rms)

# 对于长短 逻辑回归
fit <- rms::lrm(
  pfs ~ age + tp53,
  data = final_data
)

# 对于pfs可 线性回归或者cox回归


# all meta information ----------------------------------------------------

ids <- unique(metadata$PatientID)
all_meta <- readxl::read_excel('D:/cn1003/GT0918-CN-1003_data_listing.xlsx',
                               sheet = '人口统计学'
                               )
foo <- all_meta %>% filter(`Subject name or identifier` %in% ids) %>% 
  mutate(`年龄(系统产生)` = as.numeric(`年龄(系统产生)`))

write_delim(foo, 'D:/cn1003/meibo.txt', delim = '\t')




# AR Alteration survival  plot --------------------------------------------
# 还是以入组患者为分析出发点
ar_cnv_id <- (df2 %>% filter(StudyVisit=='C1D1') %>% 
  filter(Hugo_Symbol=='AR') %>% 
  select(SubjectID) %>% distinct())$SubjectID

ar_mutation_id <- (df %>% filter(StudyVisit=='C1D1') %>% 
  filter(GeneSymbol=='AR') %>% 
    filter(VariantFreq>0) %>% 
    select(SubjectID) %>% distinct())$SubjectID

ar_v7_id <- (df3 %>% filter(StudyVisit=='C1D1') %>% 
  filter(NameAST %in% c('AR-V7')) %>% select(SubjectID) %>% distinct())$SubjectID

cnv_id <- unique(c(ar_cnv_id,ar_mutation_id))
# cnv_id <- cnv_id[-which(cnv_id %in% c('02003','14002'))]

m1 <- metadata %>% select(PatientID,time, statues) %>% 
  # filter(!PatientID %in% first_id) %>%
  distinct() %>% 
  mutate(gene=ifelse((PatientID %in% cnv_id),'mutant', 'wildtype'))

fit <- survfit(Surv(time, statues) ~ gene, data = m1)
ggsurvplot(fit, pval = T)


# 数数
ar_cnv_id <- (df2 %>% 
                filter(Hugo_Symbol=='TP53') %>% 
                select(SubjectID) %>% distinct())$SubjectID

ar_mutation_id <- (df  %>% 
                     filter(GeneSymbol=='TP53') %>% 
                     filter(VariantFreq>0) %>% 
                     select(SubjectID) %>% distinct())$SubjectID


# 按照不同分组进行展开 --------------------------------------------------------------

id_b <- (metadata %>% filter(drug=='B') %>% 
  select(PatientID) %>% distinct())$PatientID

id_c <- (metadata %>% filter(drug=='C') %>% 
           select(PatientID) %>% distinct())$PatientID

id_bc <- (metadata %>% filter(drug=='C+B') %>% 
           select(PatientID) %>% distinct())$PatientID

# 批次ID
first_190329 <- readxl::read_excel('D:/cn1003/Kintor-GT0918-CN-1003_Data_03292019.xlsx',
                                   sheet = 'cfDNA_QC'
                                   )
first_id <- first_190329 %>% 
  separate(SampleID, into = c('PatientID','StudyVisit')) %>% 
  select(PatientID) %>% distinct() %>% pull(PatientID)


filter_gene <- function(gene, visit='C1D1'){
  res <- df %>% filter(GeneSymbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, VariantFreq) %>% 
    select(SubjectID,VariantFreq) %>% 
    group_by(SubjectID) %>% 
    slice_max(VariantFreq) %>% 
    dplyr::rename(!!sym(gene) := VariantFreq)
  res
}
mutation <- filter_gene('TP53')
mutation_ar <- filter_gene('AR')
mutation_atm <- filter_gene('ATM')
mutation_pik3ca <- filter_gene('PIK3CA')

ar_cnv <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')


final_data <- metadata %>% select(PatientID,time) %>% 
  left_join(mutation, by = c('PatientID'='SubjectID')) %>% 
  left_join(mutation_atm, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_pik3ca, by = c('PatientID'='SubjectID')) %>%
  left_join(mutation_ar, by = c('PatientID'='SubjectID')) %>%
  left_join(ar_cnv, by = c('PatientID'='SubjectID')) %>% 
  filter(PatientID %in% id_bc)

M <- final_data %>% 
  select(-time) %>%
  column_to_rownames('PatientID')

ha <- columnAnnotation(bar = anno_barplot(final_data$time))

Heatmap(t(M),
        name = 'ratio',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        # col = circlize::colorRamp2(c(0, 2, 4), c("navy", "white", "firebrick3")),
        col = circlize::colorRamp2(c(0,10), c("white", "firebrick3")),
        width = ncol(M)*unit(20, "mm"),
        height = nrow(M)*unit(3, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)







# mutation $ CNV ----------------------------------------------------------

# 统计同时发生CNV 和 mutation的样本

find_overlap <- function(gene, study='C1D1', ...){
  cnv_id <- df2 %>% 
    filter(StudyVisit==study) %>% 
    filter(Hugo_Symbol==gene) %>% 
    pull(SubjectID)
  
  df %>% filter(SubjectID %in% cnv_id) %>% 
    filter(GeneSymbol == gene) %>% 
    pull(SubjectID) %>% unique()
}

find_overlap('AR')




