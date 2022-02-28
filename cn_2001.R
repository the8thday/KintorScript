# for CN2001

library(survival)
library(survminer)
library(tidyverse)
library(paletteer)
library(magrittr)
library(ComplexHeatmap)


# read files --------------------------------------------------------------

# cfdna_snv1 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_data_10132019.xlsx',
#                                 sheet = 'cfDNA_SNV')
# cfdna_snv2 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_Data_08082019.xlsx',
#                                  sheet = 'cfDNA_SNV')
# cfdna_snv <- bind_rows(cfdna_snv1, cfdna_snv2)
# 
# cfdna_cnv <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_data_10132019.xlsx',
#                                 sheet = 'cfDNA_CNV')
# cfdna_cnv2 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_Data_08082019.xlsx',
#                                  sheet = 'cfDNA_CNV')
# cfdna_cnv %<>%  bind_rows(cfdna_cnv2)
# 
# cfrna_splicing <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_data_10132019.xlsx',
#                                 sheet = 'cfRNA_Splicing')
# cfrna_splicing2 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_Data_08082019.xlsx',
#                                       sheet = 'cfRNA_Splicing')
# cfrna_splicing %<>% bind_rows(cfrna_splicing2) %>% 
#   unite('ExternalID', SubjectID:StudyVisit, remove = F)
# 
# cfdna_qc1 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_data_10132019.xlsx',
#                                sheet = 'cfDNA_Summary') %>% 
#   select(-Proceeded) %>% 
#   mutate(cfDNAYield=as.double(cfDNAYield))
# cfdna_qc2 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_Data_08082019.xlsx',
#                                 sheet = 'cfDNA_Summary') %>% 
#   unite('ExternalID', c(SubjectID,StudyVisit), remove = F)
# cfdna_qc <- bind_rows(cfdna_qc1, cfdna_qc2)
# visdat::vis_dat(cfdna_qc)
# glimpse(cfdna_qc)

cfdna_snv <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_cfDNA_cfRNA_cumulative_data_20210609.xlsx',
                                sheet = 'DNA_SNV')
cfdna_cnv <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_cfDNA_cfRNA_cumulative_data_20210609.xlsx',
                                sheet = 'DNA_CNV')
cfdna_qc <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_cfDNA_cfRNA_cumulative_data_20210609.xlsx',
                                sheet = 'DNA_QC_Summary')
cfrna_splicing <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_cfDNA_cfRNA_cumulative_data_20210609.xlsx',
                                sheet = 'RNA_Splicing')
cfrna_expre <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_cfDNA_cfRNA_cumulative_data_20210609.xlsx',
                                sheet = 'RNA_Expression')


metadata <- readxl::read_excel('D:/2001/cn2001_metadata.xlsx',
                               sheet = 'Sheet1') %>% 
  arrange(days)

df1 <- cfdna_snv %>% left_join(metadata, by = c('SubjectID'='SubjectID'))
df1 %<>% unite(col = 'ExternalID', c(SubjectID, StudyVisit), remove = F)

df1 %<>% mutate(time_d = if_else(StudyVisit=='Screening'|StudyVisit=='C1D1', days,
                                 if_else(StudyVisit=='C2D1', days-28,
                                         if_else(StudyVisit=='C2D28'|StudyVisit=='C3D1', days-28*2, 
                                                 if_else(StudyVisit=='C4D27', days-28*3-27,
                                                         if_else(StudyVisit=='C4D28'|StudyVisit=='C5D1',
                                                                 days-28*4, 
                                                                 if_else(StudyVisit=='C6D28',
                                                                         days-28*6,
                                                                         if_else(StudyVisit=='EOT',
                                                                                 0,0)))
                                                 )))
))  %>% arrange(time_d) %>% 
  mutate(time_d=ifelse(time_d<0,0,time_d))

# QC
table(cfdna_qc$LibraryQC)
table(cfdna_qc$NGSQC)
cfdna_qc %>% filter(NGSQC %in% c('FAIL','AT RISK'))
cfdna_qc %>% filter(LibraryQC %in% c('FAIL','AT RISK'))

# 用药时长统计
metadata %<>% filter(days>100)
Hmisc::describe(metadata$days)
describer::describe(metadata$days)
ggplot(metadata, aes(y=days)) + geom_boxplot() + theme_bw()
ggplot(metadata, aes(x=days)) + geom_histogram() + theme_bw() +
  geom_density() +
  theme(axis.text = element_text(size = 16))

x1 <- Rmisc::CI(metadata$days,
          ci = 0.95
          )

x <- metadata$days
dat<-with(density(x),data.frame(x,y))
dat1<-dat[dat$x>x1[3]&dat$x<x1[1],]

ggplot(metadata,aes(days))+
  geom_density(fill="grey")+
  geom_vline(xintercept = x1[1],lty="dashed")+
  geom_vline(xintercept = x1[3],lty="dashed")+
  geom_area(data=dat1,aes(x=x,y=y),fill="red")+
  geom_vline(xintercept = x1[2],lty="dashed")+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,0.02))+
  theme_minimal() +
  annotate('text',x=250, y=0.016, label='[]')


# overall mutation, 统计各个mutation的个数分布
a1 <- df1 %>% 
  # filter(StudyVisit=='C1D1') %>% 
  filter(SubjectID %in% f60_id) %>% 
  unite(col = 'mut', c(Hugo_Symbol,HGVSp), remove = F) %>% 
  group_by(mut) %>% count() %>% arrange(desc(n))
# 患者在不同时期测序的情况
a1 <- cfdna_qc %>% dplyr::select(SubjectID, StudyVisit, LibraryYield) %>% 
  filter(SubjectID %in% f60_id) %>% 
  pivot_wider(names_from = StudyVisit, values_from = LibraryYield) %>% 
  mutate(across(where(is.double), ~if_else(is.na(.), 'NA', 'YES')))
write_delim(a1, 'D:/2001/foo.txt', delim = '\t')


# survival analysis -------------------------------------------------------

fit <- survfit(Surv(days, status) ~ dosing, data = metadata)
summary(fit)
surv_summary(fit)
ggsurvplot(fit, pval = T)
# log rank value
survdiff1 <- survdiff(Surv(days, status) ~ dosing, data = metadata)
plot1 <- ggsurvplot(fit, pval = T)
plot1$plot

fit2 <- survfit(Surv(days, status) ~ subtype, data = metadata)
ggsurvplot(fit2, pval = T)

# cox analysis
reg <- coxph(Surv(days, status) ~ AR_percent, data = metadata)
summary(reg)
survminer::ggforest(reg, data = metadata)

df_reg <- metadata
df_reg$dosing  <- factor(df_reg$dosing, levels = c(200,300), 
                         labels = c('200mg', '300mg'))
reg2 <- coxph(Surv(days, status) ~ AR_percent + subtype + dosing, 
              data = df_reg)
survminer::ggforest(reg2, data = df_reg)
cox.zph(reg2)

# AR percent
mm <- metadata %>% mutate(arP=if_else(AR_percent>0.7,'ARH','ARL'))
mm <- mm %>% filter(SubjectID %in% f60_id)
fit <- surv_fit(Surv(days, status)  ~ arP, data = mm)
ggsurvplot(fit, pval = T)



# QC Analysis -------------------------------------------------------------
library(patchwork)
library(ggbreak)

qc1 <- cfdna_qc %>% left_join(metadata, by = 'SubjectID') %>% 
  arrange(days)
qc1.c1 <- qc1 %>% filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  mutate(tmb=case_when(
    TMBScore < 14 ~ 'low',
    TMBScore > 14 ~ 'high',
    TRUE ~ 'NA'
  )) %>% 
  filter(!is.na(TMBScore))
qc1.c1 %<>% filter(SubjectID %in% hr_id)

fit <- surv_fit(Surv(days, status)~tmb, data = qc1.c1)
ggsurvplot(fit,
           conf.int = F,
           pval = T,
           risk.table = T,
           surv.median.line = "hv",
           pval.method = T
           )

reg <- coxph(Surv(days, status) ~ TMBScore, data = qc1.c1)
ggforest(reg,
         data = qc1.c1,
         main = "Hazard ratio"
         )

# 对于DNA Yield
zz <- qc1 %>% select(SubjectID, StudyVisit, DNAYield, days) %>% 
  arrange(StudyVisit) %>% 
  pivot_wider(names_from = StudyVisit, values_from = DNAYield) %>% 
  arrange(days) %>% 
  mutate(C1D1=if_else(is.na(C1D1),Screening,C1D1)) %>% 
  mutate(C2D28=if_else(is.na(C2D28),C3D1,C2D28)) %>% 
  mutate(C5D1=if_else(is.na(C5D1),C4D28,C5D1)) %>% 
  mutate(C5D1=if_else(is.na(C5D1),C4D27, C5D1)) %>% 
  mutate(EOT=if_else(
    is.na(EOT),if_else(is.na(C6D28),if_else(is.na(C5D1),
                                            if_else(is.na(C2D28),if_else(
                                              is.na(C2D1),C1D1,C2D1
                                            ),C2D28),C5D1),C6D28),EOT
  )) %>% 
  select(SubjectID,days,C1D1,C2D1,C2D28,C5D1,EOT)

write_delim(zz, 'D:/2001/zz.txt', delim = '\t')

zzz <- zz %>% mutate(ratio=(as.double(EOT)-as.double(C1D1))/as.double(C1D1))
zzz$SubjectID <- factor(zzz$SubjectID, levels = zzz$SubjectID)
zzz %<>% mutate(across(starts_with(c('E','C')), as.double)) 
  # filter(SubjectID !='B02007')
zzz %<>% filter(SubjectID %in% f60_id)

p1 <- ggplot(zzz, aes(x=SubjectID,y=C1D1)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size = 8))
p2 <- ggplot(zzz, aes(x=SubjectID,y=EOT)) +
  geom_col() +
  # scale_y_break(c(200, 500), scales = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size = 8))

p3 <- ggplot(zzz, aes(x=SubjectID,y=ratio)) +
  geom_col() +
  # scale_y_break(c(1, 5), scales = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size = 8))

p1/p2/p3

library(aplot)
ap <- p3 %>% 
  insert_top(p2) %>% 
  insert_top(p1)

# 生存期较长的几个样本
describer::describe(metadata$days)
m200 <- metadata %>% filter(days>200)



# 统计各基因突变 -----------------------------------------------------------------
df1 %>% filter(Hugo_Symbol=='TP53') %>% 
  group_by(SubjectID) %>% count() %>% 
  arrange(desc(n))

gene_mutation_count <- function(gene){
  df1 %>% filter(Hugo_Symbol==gene) %>% 
    select(HGVSp,ExternalID, Clinvar,days, time_d, status,
           AltDepth, Variant_Value, Mutation_Status) %>% 
    # filter(Variant_Value>1) %>%
    arrange(days) %>% 
    distinct()
}
bar <- gene_mutation_count('PIK3CA') %>% data.frame()

write_delim(bar, 'D:/2001/mutation.txt', delim = '\t')


# C1D1期检测到某基因出现的数目
df1 %>% filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  filter(Hugo_Symbol=='ERBB2') %>% 
  group_by(SubjectID) %>% 
  count()
# 统计某基因所有突变 在 不同采样期出现的数目
aaaa <- df1 %>% filter(Hugo_Symbol=='') %>% 
  filter(!SubjectID %in% f60_id) %>% 
  select(SubjectID, ExternalID, StudyVisit,Variant_Value, HGVSp) %>% 
  distinct() %>% group_by(HGVSp, StudyVisit) %>% 
  summarise(n = n()) %>% arrange(n) %>% 
  pivot_wider(names_from = StudyVisit, 
              values_from = n
              )
write_delim(aaaa, 'D:/2001/genemutation.txt', delim = '\t')

# 野生型和突变型 生存曲线(C1D1期)
c1d1.samples <- df1 %>% filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp, days, status) %>% 
  filter(Hugo_Symbol=='ESR1') %>% select(SubjectID) %>% distinct()

m1 <- df1 %>% select(SubjectID, days, status) %>% distinct() %>%
  filter(SubjectID %in% f60_id) %>% 
  mutate(ESR1=if_else(SubjectID %in% c1d1.samples$SubjectID, 'mutant',
                      'wildtype'))
fit <- survfit(Surv(days, status) ~ ESR1, data = m1)
ggsurvplot(fit, pval = T)

# 只考虑C1D1期 pathogenic
c1d1.samples <- df1 %>% 
  filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  filter(Clinvar %in% c('Pathogenic', 'Likely_pathogenic', 
                        'Pathogenic/Likely_pathogenic')) %>% 
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp, days, status) %>% 
  filter(Hugo_Symbol=='TP53') %>% select(SubjectID) %>% distinct()



# plot --------------------------------------------------------------------

aa <- metadata %>% select(SubjectID, days, subtype) %>% 
  distinct() %>% arrange(days)

aa <- metadata %>% 
  filter(SubjectID %in% f60_id)
aa$SubjectID <- factor(aa$SubjectID, levels = aa$SubjectID)

ggplot(aa) + geom_bar(aes(x = SubjectID, y = days, fill=subtype), 
                            stat = 'identity') + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  ggsci::scale_fill_aaas()


# C1D1期某个基因的突变频率，突变频率最高的mutation
c('C1D1','Screening')
c('C2D28', 'C3D1')
eot <- c('EOT','C6D28','C5D1','C4D27', 'C4D28')
filter_gene <- function(gene, visit=c('C1D1','Screening')){
  res <- df1 %>% filter(Hugo_Symbol %in% gene) %>% 
    filter(StudyVisit %in% visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}

mut_TP53 <- filter_gene('TP53')
mut_PIK3CA <- filter_gene('PIK3CA')
mut_ESR1 <- filter_gene('ESR1')
mut_AKT1 <- filter_gene('AKT1')
mut_CHEK2 <- filter_gene('CHEK2')
mut_ATM <- filter_gene('ATM')

FGFR1_cnv <- df2 %>% filter(Hugo_Symbol=='FGFR1') %>% 
  filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('FGFR1_CNV'='Variant_Value')

final_data <- df1 %>% select(SubjectID,days) %>% 
  distinct() %>% 
  arrange(days) %>% 
  left_join(mut_TP53, by = c('SubjectID')) %>% 
  left_join(mut_PIK3CA, by = c('SubjectID')) %>%
  left_join(mut_ESR1, by = c('SubjectID')) %>%
  left_join(mut_AKT1, by = c('SubjectID')) %>%
  left_join(mut_CHEK2, by = c('SubjectID')) %>% 
  left_join(mut_ATM, by = c('SubjectID')) %>% 
  left_join(FGFR1_cnv, by = c('SubjectID'))

final_data %<>% filter(SubjectID %in% f60_id)

M <- final_data %>% 
  select(-days) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$days))

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
        col = circlize::colorRamp2(c(0,5), c("white", "firebrick3")),
        width = ncol(M)*unit(20, "mm"),
        height = nrow(M)*unit(6, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = F,
        # column_names_gp = gpar(fontsize=1)
)

# 某基因的全部突变

filter_aa <- function(gene, aa){
  res <- df1 %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit %in% c('C1D1','Screening')) %>%
    filter(HGVSp==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}
filter_aa('PIK3CA','p.His1047Arg')


all_gene <- c(
  "p.Val812Ile",
  # "p.Arg1018His",
  # "p.Pro1140Ala",
  # "p.Thr729Ala",
  "p.Gly697Ala",
  "p.Asp739Tyr",
  "p.Leu725Ser"
)

final_data <- metadata %>% select(SubjectID,days) %>% 
  arrange(days)
for(i in all_gene){
  foo <- filter_aa('PIK3CA', i)
  final_data <- final_data %>%
    left_join(foo, by = c('SubjectID'='SubjectID'))
  
}
final_data %<>% filter(SubjectID %in% f60_id)
M <- final_data %>% 
  select(-days) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$days))


# CNV   -------------------------------------------------------------

sort(table(cfdna_cnv$Hugo_Symbol)) # FGFR1,CCND1,AKT3,MYC,NF2,CDKN2A
df2 <- cfdna_cnv %>% left_join(metadata, by = c('SubjectID'))

# AR CNV
cfdna_cnv %>% filter(Hugo_Symbol=='AR') #并没有AR CNV
cfdna_cnv %>% filter(Hugo_Symbol=='CDKN2A')


# 分时期对某些基因CNV 进行展示
i_gene <- Hmisc::Cs(FGFR1,CCND1,AKT3,MYC,NF2,CDKN2A)
i_cnv <- df2 %>% filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(days) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))
# i_cnv <- i_cnv[, order(names(i_cnv))]
ff <- metadata %>% select(SubjectID,days) %>% 
  arrange(days) %>% 
  left_join(i_cnv, by = c('SubjectID')) %>% 
  filter(SubjectID %in% f60_id)
M <- ff %>% 
  select(-days) %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(bar = anno_barplot(ff$days))


# C1D1期AR CNV 出现与否的survival
cnv.s <- function(gene){
  c1d1 <- df2 %>% filter(StudyVisit %in% c('C1D1','Screening')) %>% 
    filter(Hugo_Symbol==gene) %>% 
    select(SubjectID) %>% distinct()
  c1d1$SubjectID
}

c1d1.FGFR1 <- cnv.s('FGFR1')
c1d1.CCND1 <- cnv.s('CCND1')
c1d1.AKT3 <- cnv.s('AKT3')
c1d1.MYC <- cnv.s('MYC')
c1d1.NF2 <- cnv.s('NF2')

m1 <- df1 %>% select(SubjectID, days, status) %>% distinct() %>%
  mutate(FGFR1=if_else(SubjectID %in% c1d1.FGFR1, 2,1)) %>% 
  mutate(CCND1=if_else(SubjectID %in% c1d1.CCND1, 2,1)) %>% 
  mutate(AKT3=if_else(SubjectID %in% c1d1.AKT3, 2,1)) %>% 
  mutate(MYC=if_else(SubjectID %in% c1d1.MYC, 2,1)) %>% 
  mutate(NF2=if_else(SubjectID %in% c1d1.NF2, 2,1))

reg <- coxph(Surv(days, status) ~ FGFR1+CCND1+AKT3+MYC+NF2,
             data = m1)
summary(reg)
survminer::ggforest(reg, data = m1)

# logitic regression ------------
m1 %<>% mutate(class=if_else(days>150, 1, 0)) 
m2 <- m1 %>% select(-SubjectID, -days, -status)

lr_model <- glm(class~., data = m2, family = binomial)
lm_model <- lm(formula = days ~ FGFR1+CCND1+AKT3+MYC+NF2, data = m1)
plot(performance::check_collinearity(lr_model))
summary(lr_model)
summary(lm_model)




# splicing ----------------------------------------------------------------
sort(table(cfrna_splicing$NameAST))
# 发生频率最高的为AR-2_3、 AR-4_5、AR-3_4、AR-FL-7_8

df3 <- cfrna_splicing %>% left_join(metadata, 
                                    by = c('SubjectID')) %>% 
  unite(ExternalID, c(SubjectID,StudyVisit), remove = FALSE)

arv7 <- df3 %>% 
  filter(NameAST %in% c('AR-V7',
                        'AR-V3', 'AR-V9'
  )) %>% 
  filter(SupportReads >= 3) %>% 
  select(SubjectID, ExternalID, StudyVisit, 
         NameAST, SupportReads, days, status) %>% 
  arrange(desc(days))

write_delim(arv7, 'D:/2001/foo.txt', delim = '\t')

# C1D1
df3 %>% filter(StudyVisit %in% c('C1D1','Screening'))
df3 %>% filter(NameAST=='AR-FL-7_8') %>% 
  filter(SubjectID %in% c('B02007'))


# Expression --------------------------------------------------------------

m1 <- metadata %>% dplyr::select(SubjectID, days)
cfrna_expre1 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_data_10132019.xlsx',
                                     sheet = 'cfRNA_Expression')
cfrna_expre2 <- readxl::read_excel('D:/2001/Kintor-GT0918-CN-2001_Data_08082019.xlsx',
                                      sheet = 'cfRNA_Expression')
# cfrna_expre <- bind_rows(cfrna_expre1, cfrna_expre2) %>% 
#   left_join(m1, 
#             by = c('SubjectID'))
delet.genes <- setdiff(unique(names(cfrna_expre2)), unique(names(cfrna_expre1)))
cfrna_expre2_id <- cfrna_expre2$RequisitionID %>% 
  str_c('_R')

cfrna_expre %<>% left_join(m1, by = c('SubjectID')) %>% 
  dplyr::select(-{{ delet.genes }})
# cfrna_expre %>% select(where(~sum(!is.na(.x))>20))
cfrna_expre <- cfrna_expre[, !apply(cfrna_expre,2, function(x){sum(is.na(x))>20})] %>% 
  mutate(pici = if_else(RequisitionID %in% cfrna_expre2_id, 'a', 'b'))

#-----------去除批次效应
library(limma)
cfrna_num <- cfrna_expre %>% 
  select(-pici, -c(StudyVisit,SubjectID,days)) %>% 
  column_to_rownames('RequisitionID') %>% t() %>% 
  as.data.frame()

design <- model.matrix(~ 0 + pici)

df_rBE <- removeBatchEffect(
  cfrna_num,
  batch = cfrna_expre$pici
)

mmm <- cfrna_expre %>% select(StudyVisit,SubjectID,RequisitionID,days,pici)

cfrna_expre <- data.frame(t(df_rBE)) %>% rownames_to_column() %>% 
  as_tibble() %>% 
  left_join(mmm, by = c('rowname'='RequisitionID')) %>% 
  rename(RequisitionID=rowname)


# select <- dplyr::select
# 不同时期的基因表达
fpkm_c1d1 <- cfrna_expre %>% 
  filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  select(-c(StudyVisit,RequisitionID,days,pici)) %>% 
  column_to_rownames('SubjectID') %>% t() %>% 
  as.data.frame()
visdat::vis_miss(fpkm_c1d1)
fpkm_c1d1 %<>% replace(is.na(.), 0)

col_anno <- cfrna_expre %>% 
  filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  select(SubjectID, days) %>% 
  distinct() %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(time = anno_barplot(col_anno))

Heatmap(
  log10(fpkm_c1d1 + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=12),
  row_names_gp = gpar(fontsize=0.5),
  cluster_rows = T,
  cluster_columns = T
)

# 层次聚类
hc <- hclust(dist(log10(fpkm_c1d1 + 1),method = "euclidean"),
             method = "complete")
summary(hc)
plot(hc, labels = F
     )
# heatmap(as.matrix(hc))
rect.hclust(hc, k=5)
genes <- as.data.frame(cutree(hc, k=5))
names(genes) <- 'class'
gg <- genes %>% filter(class==2) %>% 
  rownames_to_column()
write_delim(gg, 'D:/2001/hc.txt', delim = '\t')


dend1 <- as.dendrogram(hc)
plot(dend1,
     edgePar = list(lab.cex = 0.1),
     nodePar = list(lab.cex = 0.5)
)
# 特定基因的表达差异
i_gene <- Hmisc::Cs(AR,TP53,PIK3CA,ESR1,CHEK2,AKT1,ERBB2,MYC)
fpkm_c1d1.i <- fpkm_c1d1[i_gene, ]
Heatmap(
  log10(fpkm_c1d1.i + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=12),
  row_names_gp = gpar(fontsize=12),
  cluster_rows = T,
  cluster_columns = T
)




# 基因表达值 做线性回归 预测


# AR GENE & AR percent
ar_p <- cfrna_expre %>% left_join(metadata, by = c('SubjectID')) %>% 
  filter(StudyVisit %in% c('C1D1','Screening')) %>% 
  select(SubjectID, AR, AR_percent, days.x) %>% 
  arrange(days.x)


arpp <- ar_p %>% pivot_longer(cols = -SubjectID)
arpp$SubjectID <- factor(arpp$SubjectID, levels = ar_p$SubjectID)

arpp %>% 
  ggplot(., aes(x=SubjectID, y=value)) + 
  geom_bar(stat = 'identity') +
  facet_grid(name~., scales = 'free') +
  ggthemes::theme_calc() +
  theme(axis.text.x = element_text(
    angle = 90
  ))

cc <- correlation::cor_test(ar_p, 
                      x = 'AR',
                      y = 'AR_percent',
                      method = 'spearman'
                      )

cc %>% 
  summary(redundant=T) %>% 
  plot()

write_delim(as_tibble(cc), 
            file = 'D:/2001/foo.txt',
            delim = '\t'
            )

correlation::cor_test(log(ar_p$AR), ar_p$AR_percent)
cor.test(log(ar_p$AR), ar_p$AR_percent)


## AR 表达随时间的变化
fpkm <- cfrna_expre %>% 
  select(-c(SubjectID,StudyVisit,days)) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames('RequisitionID') %>% 
  t() %>% as.data.frame()

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

tpms <- apply(fpkm,2,fpkmToTpm) %>% 
  as.data.frame() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
kk <- cfrna_expre %>% select(RequisitionID,
                             SubjectID,
                             StudyVisit,
                             days
                             )
cfrna_expre <- tpms %>% left_join(kk,
  by = c('rowname'='RequisitionID')
)

ar_expre <- cfrna_expre %>% 
  arrange(days, StudyVisit) %>% 
  select(SubjectID, StudyVisit, days, AR) %>% 
  distinct() %>% 
  pivot_wider(names_from = StudyVisit, values_from = AR) %>% 
  mutate(C1D1=if_else(is.na(C1D1),Screening,C1D1)) %>% 
  mutate(C2D28=if_else(is.na(C2D28),C3D1,C2D28)) %>% 
  mutate(C5D1=if_else(is.na(C5D1),C4D28,C5D1)) %>% 
  mutate(C5D1=if_else(is.na(C5D1),C4D27, C5D1)) %>% 
  mutate(EOT=if_else(
    is.na(EOT),if_else(is.na(C6D28),if_else(is.na(C5D1),
                                            if_else(is.na(C2D28),if_else(
                                              is.na(C2D1),C1D1,C2D1
                                            ),C2D28),C5D1),C6D28),EOT
  )) %>% 
  select(SubjectID,days,C1D1,C2D1,C2D28,C5D1,EOT)

M <- ar_expre %>% 
  select(-days) %>% 
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(ar_expre$days))

Heatmap(log10(t(M)+1),
        name = 'log10(TPM)',
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
        width = ncol(M)*unit(40, "mm"),
        height = nrow(M)*unit(2, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        column_names_gp = gpar(fontsize=8)
)


# gene survival -----------------------------------------------------------
aa <- Hmisc::Cs(
  p.Phe134Cys,
  p.Tyr220Cys,
  p.Arg306Ter,
  p.Pro177_Cys182del,
  p.Leu194His,
  p.Arg267Trp,
  p.Leu265del,
  p.Lys351Ter,
  p.Arg175His,
  p.Val143Met,
  p.Gln144Pro,
  p.Arg213Leu,
  p.Asp281Gly,
  p.Lys132Arg,
  p.Glu358Val
)
aa <- Hmisc::Cs(
  p.Met983Ile,
  p.Glu545Lys,
  p.Tyr343His,
  p.Gly118Asp,
  p.Lys711Thr,
  p.Glu453Lys,
  p.Asn345Lys
)
# 按照SubjectID
df_p53 <- df1 %>% mutate(stat_p53 = ifelse(GeneSymbol=='PIK3CA' & 
                                            HGVSp %in% aa,
                                          'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)

d_tp53 <- metadata %>% right_join(df_p53, by = c('SubjectID'='SubjectID'))

fit <- survfit(Surv(days, status) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)


# 按照ExternalID
foo <- df1 %>% select(ExternalID, days, time_d, status) %>% distinct()
df_p53 <- df1 %>% mutate(stat_p53 = ifelse(GeneSymbol=='TP53' & 
                                            HGVSp %in% aa,
                                          'yes', 'no')) %>%
  # select(SubjectID, stat_p53, ExternalID, time, time_d, stat) %>%
  select(ExternalID, stat_p53) %>%
  distinct() %>% 
  arrange(ExternalID, desc(stat_p53)) %>% 
  group_by(ExternalID) %>% 
  top_n(1) %>% ungroup()
d_tp53 <- df_p53 %>% left_join(foo, by = c('ExternalID'))
fit <- survfit(Surv(time_d, status) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

write_delim(d_tp53, 'D:/2001/d_tp53.txt', delim = '\t')
d_tp53 <- read_delim('D:/2001/first_p53.txt', delim = '\t')


# 只考虑患者入组时的状态
ruzu <- df1 %>% filter(StudyVisit=='C1D1')
df_p53 <- ruzu %>% mutate(stat_p53 = ifelse(GeneSymbol=='PIK3CA' & 
                                              HGVSp %in% aa,
                                            'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)

d_tp53 <- metadata %>% right_join(df_p53, by = c('SubjectID'='SubjectID'))
fit <- survfit(Surv(days, status) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)


######################### 从ARV7 出现与否的角度进行分析#########################
# 首先是按照ExternalID进行
foo <- df3 %>% select(ExternalID, time, time_d, stat) %>% distinct()
df_arv7 <- df3 %>% mutate(stat_arv7 = if_else(NameAST=='AR-V7', 'positive', 'negtive')) %>% 
  select(ExternalID, stat_arv7) %>%
  distinct() %>% 
  arrange(ExternalID, desc(stat_arv7)) %>% 
  group_by(ExternalID) %>% 
  top_n(1) %>% ungroup()
d_arv7 <- df_arv7 %>% left_join(foo, by = c('ExternalID'))
d_arv7 %<>% filter(time_d > 0)
fit <- surv_fit(Surv(time_d, stat) ~ stat_arv7, data = d_arv7)
ggsurvplot(fit, pval = T)



























# 去掉HER+++ ----------------------------------------------------------------

her_id <- (metadata %>% filter(subtype == 'HER+++'))$SubjectID
hr_id <- (metadata %>% filter(subtype == 'HR+'))$SubjectID
trineg_id <- (metadata %>% filter(subtype == '三阴性'))$SubjectID

m_her <- metadata %>% filter(subtype != 'HER+++')
  
  
fit2 <- survfit(Surv(days, status) ~ subtype, data = m_her)
ggsurvplot(fit2, pval = T)

# 去掉60天以内样本的分析
f60_id <- (metadata %>% filter(days > 100))$SubjectID



#
# 重新定义EOT -----------------------------------------------------------------
a1 <- cfdna_qc %>% dplyr::select(SubjectID, StudyVisit, DNAYield) %>% 
  pivot_wider(names_from = StudyVisit, values_from = DNAYield) %>% 
  mutate(across(where(is.double), ~if_else(is.na(.), 'NA', 'YES')))

eot <- c('EOT','C6D28','C5D1','C4D27')
eot_c4d28_id <- c('B05009','B01010','B01007')
eot_c2d28_id <- c('B04007','B01004','B01005','B05012','B01014','B01015','B01003')
eot_c2d1_id <- c('B02007')

filter_gene <- function(gene, visit=eot){
  res <- df1 %>% filter(Hugo_Symbol %in% gene) %>% 
    filter(StudyVisit %in% visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}


filter_gene_id <- function(gene, id, visit=eot){
  res <- df1 %>% filter(Hugo_Symbol %in% gene) %>% 
    filter(StudyVisit %in% visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}

ff <- function(gene){
  foo <- filter_gene(gene) %>% 
  bind_rows(filter_gene_id(gene, eot_c4d28_id, c('C4D28'))) %>% 
  bind_rows(filter_gene_id(gene, eot_c2d28_id, c('C2D28'))) %>% 
  bind_rows(filter_gene_id(gene, eot_c2d1_id, c('C2D1')))
  foo
  }


mut_TP53 <- ff('TP53')
mut_PIK3CA <- ff('PIK3CA')
mut_ESR1 <- ff('ESR1')
mut_AKT1 <- ff('AKT1')
mut_CHEK2 <- ff('CHEK2')
mut_ATM <- ff('ATM')

FGFR1_cnv <- df2 %>% filter(Hugo_Symbol=='FGFR1') %>% 
  filter(StudyVisit %in% eot) %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('FGFR1_CNV'='Variant_Value')

filter_cnv <- function(id, visit, gene='FGFR1'){
  df2 %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit %in% visit) %>% 
    select(SubjectID, Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    rename('FGFR1_CNV'='Variant_Value')
}

FGFR1_cnv %<>% 
  bind_rows(filter_cnv(eot_c4d28_id, c('C4D28'))) %>% 
  bind_rows(filter_cnv(eot_c2d28_id, c('C2D28'))) %>% 
  bind_rows(filter_cnv(eot_c2d1_id, c('C2D1')))

final_data <- df1 %>% select(SubjectID,days) %>% 
  distinct() %>% 
  arrange(days) %>% 
  left_join(mut_TP53, by = c('SubjectID')) %>% 
  left_join(mut_PIK3CA, by = c('SubjectID')) %>%
  left_join(mut_ESR1, by = c('SubjectID')) %>%
  left_join(mut_AKT1, by = c('SubjectID')) %>%
  left_join(mut_CHEK2, by = c('SubjectID')) %>% 
  left_join(mut_ATM, by = c('SubjectID')) %>% 
  left_join(FGFR1_cnv, by = c('SubjectID'))

final_data %<>% filter(SubjectID %in% f60_id)

M <- final_data %>% 
  select(-days) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$days))

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
        col = circlize::colorRamp2(c(0,5), c("white", "firebrick3")),
        width = ncol(M)*unit(20, "mm"),
        height = nrow(M)*unit(6, "mm"),
        # top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)


# 某基因的全部突变
filter_aa <- function(gene, aa){
  res <- df1 %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit %in% eot) %>%
    filter(HGVSp==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}

filter_aa_id <- function(gene, aa, id, visit){
  res <- df1 %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit %in% visit) %>%
    filter(HGVSp==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}

fff <- function(gene, aa){filter_aa(gene = gene, aa = aa) %>% 
  bind_rows(filter_aa_id(gene = gene, aa = aa, 
             id = eot_c4d28_id, visit = c('C4D28'))) %>% 
  bind_rows(filter_aa_id(gene = gene, aa = aa, 
                         id = eot_c2d28_id, visit = c('C2D28'))) %>% 
  bind_rows(filter_aa_id(gene = gene, aa = aa, 
                         id = eot_c2d1_id, visit = c('C2D1')))}

all_gene <- c(
  "p.His1047Arg",
  "p.Met983Ile",
  "p.Glu545Lys",
  "p.Tyr343His",
  "p.Gly118Asp",
  "p.Lys711Thr",
  "p.Glu453Lys",
  "p.His1047Leu",
  "p.Asn345Lys",
  "p.Glu542Lys",
  "p.Asn1044Lys",
  "p.Gly1049Arg",
  "p.Glu545Asp"
)

final_data <- metadata %>% select(SubjectID,days) %>% 
  arrange(days)
for(i in all_gene){
  foo <- fff('PIK3CA', i)
  final_data <- final_data %>%
    left_join(foo, by = c('SubjectID'='SubjectID'))
  
}
final_data %<>% filter(SubjectID %in% f60_id)

M <- final_data %>% 
  select(-days) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$days))


# CNV
i_gene <- Hmisc::Cs(FGFR1,CCND1,AKT3,MYC,NF2,CDKN2A)

i_cnv <- df2 %>% filter(StudyVisit %in% eot) %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(days) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))

i_cnv_id <- function(id, visit){
  df2 %>% filter(StudyVisit %in% visit) %>% 
    filter(Hugo_Symbol %in% i_gene) %>% 
    arrange(days) %>% 
    select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
    select(sort(names(.)))
}

i_cnv %<>% 
  # bind_rows(i_cnv_id(eot_c4d28_id, c('C4D28'))) %>%
  bind_rows(i_cnv_id(eot_c2d28_id, c('C2D28'))) %>% 
  bind_rows(i_cnv_id(eot_c2d1_id, c('C2D1')))

ff <- metadata %>% select(SubjectID,days) %>% 
  arrange(days) %>% 
  left_join(i_cnv, by = c('SubjectID')) %>% 
  filter(SubjectID %in% f60_id)

M <- ff %>% 
  select(-days) %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(bar = anno_barplot(ff$days))




