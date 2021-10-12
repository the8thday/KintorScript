#
library(tidyverse)
library(ComplexHeatmap)


# read -------------------------------------------------------------------

cfdna_snv <- readxl::read_excel('D:/us1002//GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                                sheet = 'DNA_SNV')

cfdna_cnv <- readxl::read_excel('D:/us1002/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                                sheet = 'DNA_CNV')
cfdna_qc <- readxl::read_excel('D:/us1002/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                               sheet = 'DNA_QC_Summary') %>% 
  mutate(StudyVisit=if_else(StudyVisit=='End of Treatment', 'EOT', StudyVisit))
cfrna_splicing <- readxl::read_excel('D:/us1002/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                                     sheet = 'RNA_Splicing')
cfrna_snv <- readxl::read_excel('D:/us1002/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                                sheet = 'RNA_SNV_Confirm')

metadata <- readxl::read_excel('D:/us1002/GT0918-US-1002_Subject_status_2020-1014-general_dates-TH_12Mar2021.xlsx')
d <- metadata %>% 
  mutate(time = `EOT`-`Date of 1st Dosing`) %>% 
  mutate(stat = if_else(`Treatment discontunation`%in% 
                          c('Bone or Radiographic Disease Progression')|
                          `if others, specify`%in% c('Rise in PSA',
                                                     'PI decision due to rising PSA',
                                                     'rising PSA and intolerable toxicities',
                                                     'Clinical progression'
                          ), 
                        2, 1)) %>% 
  select(c('pt. ID', 'dosing', 'Prior Treat', 'time', 'stat')) %>%
  rename(drug = `Prior Treat`)
d$time %<>% as.vector

ggplot(d, aes(y=time)) + geom_boxplot() + theme_bw() +
  theme(axis.text = element_text(size = 15))
# write_delim(d, 'D:/us1002/tidy_metadata.txt', delim = '\t')

df <- cfdna_snv %>% left_join(d, by = c('SubjectID'='pt. ID')) %>% 
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
df %>% 
  # filter(StudyVisit=='C1D1') %>% 
  unite(col = 'mut', c(Hugo_Symbol,HGVSp_Short), remove = F) %>% 
  group_by(mut) %>% count() %>% arrange(desc(n)) %>% filter(n>=2)
# 患者在不同时期测序的情况
a1 <- cfdna_qc %>% select(SubjectID, StudyVisit, LibraryYield) %>% 
  pivot_wider(names_from = StudyVisit, values_from = LibraryYield) %>% 
  mutate(across(where(is.double), ~if_else(is.na(.), 'NA', 'YES')))
a2 <- cfdna_qc %>% select(SubjectID, StudyVisit, LibraryYield) %>% 
  group_by(StudyVisit) %>% count() %>% arrange(desc(n))

#QC情况
table(cfdna_qc$LibraryQC)
table(cfdna_qc$NGSQC)

df %>% select(SubjectID, dosing, drug, time, stat) %>% 
  distinct() %>% 
  write_delim('D:/us1002/for_fenzu.txt', delim = '\t')



# survival analysis -------------------------------------------------------

fit <- survfit(Surv(time, stat) ~ drug, data = d)

pp <- ggsurvplot(fit, 
           data = m1,
           surv.median.line = "hv",
           legend.title = "AR",
           legend.labs = c("Mutant+CNV", "Wildtype"),
           legend = c(0.7,0.8),
           censor.shape="|", censor.size = 4,
           pval = TRUE,
           pval.method = TRUE,
           pval.size = 8,
           # pval.coord = c(0, 0.03),
           conf.int = FALSE,
           conf.int.style = "step",
           add.all = FALSE,
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           risk.table.fontsize = 8,
           risk.table.y.text = TRUE,
           fontsize = 8,
           font.x = c(20),
           font.y = c(20),
           # font.tickslab = c(20, "plain"),
           ncensor.plot = FALSE,
           ncensor.plot.height = 0.25,
           # cumevents.y.text = 8,
           # ylab="Cumulative survival (percentage)",
           xlab = " Time (Days)",
           palette = 'aaas',
           ggtheme = cowplot::theme_cowplot(font_size = 20)
           # ggtheme = theme_survminer() +
           #   theme(legend.text = element_text(size = 20),
           #         legend.title = element_text(size = 20),
           #         axis.text.x = element_text(size = 20),
           #         text = element_text(size = 20),
           #         plot.title = element_text(hjust = 0.5)
           #   )
)
pp$plot <- pp$plot + 
  labs(title    = "Survival curves",
       subtitle = "Based on Kaplan-Meier estimates"
       )
pp

pdf("D:/us1002/plot_43_1.pdf",
    width = 10,
    height = 10
    )
print(pp, newpage = FALSE)
dev.off()

ggsave(filename = 'D:/us1002/plot_1.png',
       plot = print(pp),
       width = 10,
       height = 10,
       dpi = 500
       )



# 设置用药期长短 ----------------------------------------------------------
bar <- df %>% select(c(SubjectID, ExternalID, StudyVisit, Variant_Value, 
                       Hugo_Symbol, HGVSc, HGVSp_Short, time, stat)) %>% 
  arrange(desc(time)) %>% 
  unite('variant', c(Hugo_Symbol, HGVSp_Short), sep = ':', remove = T)

t <- bar %>% select(SubjectID, time, stat) %>% distinct() #按照t的分布情况划分
shapiro.test(t$time)
describer::describe(t$time)
# 
long <- t$SubjectID[1:15]
short <- t$SubjectID[16:44]




#
# tp53 --------------------------------------------------------------------

# 统计long/short中出现TP53的比例
count_sample_gene <- function(cut, gene){
  df %>% filter(SubjectID %in% cut) %>% filter(Hugo_Symbol==gene) %>% 
  group_by(SubjectID) %>% count()
  }
count_sample_gene(short, 'TP53')

cout_ExternalID_gene <- function(df, gene){
  df %>% filter(Hugo_Symbol==gene) %>% group_by(ExternalID) %>% 
    summarise(n = n())
}
cout_ExternalID_gene(df, 'TP53')

# 统计long/short中tp53突变情况
gene_mutation_count <- function(cut, gene){
  df %>% filter(SubjectID %in% cut) %>% filter(Hugo_Symbol==gene) %>% 
  # group_by(HGVSp_Short) %>%
  # group_by(HGVSc) %>%
  # count()
  select(HGVSp_Short,ExternalID, Clinvar, CLIN_SIG, time, time_d, 
         stat,AltDepth, Variant_Value) %>% 
    # filter(Variant_Value>0.5) %>%
    arrange(ExternalID) %>% 
    distinct()
  }
bar <- gene_mutation_count(long, 'TP53') %>% data.frame() %>% 
  distinct() %>% 
  arrange(time)

write_delim(bar, 'D:/us1002/mutation.txt', delim = '\t')

# 统计某基因所有突变 在 不同采样期出现的数目
aaaa <- df %>% filter(Hugo_Symbol=='AR') %>% 
  select(SubjectID, ExternalID, StudyVisit,Variant_Value, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short, StudyVisit) %>% 
  summarise(n = n()) %>% arrange(n) %>% 
  pivot_wider(names_from = StudyVisit, values_from = n)
write_delim(aaaa, 'D:/us1002/genemutation.txt', delim = '\t')

# C1D1期检测到某基因出现的数目
df %>% filter(StudyVisit=='C1D1') %>% 
  filter(Hugo_Symbol=='DNMT3A') %>% 
  group_by(SubjectID) %>% 
  count()
# 
hxd <- df %>% filter(Hugo_Symbol=='AR') %>% 
  # filter(StudyVisit=='C1D1') %>%
  filter(HGVSp_Short=='p.T878S') %>% 
  select(SubjectID) %>% distinct()

# TP53 野生型和突变型 生存曲线(C1D1期)
c1d1.samples <- df %>% filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, time, stat) %>% 
  filter(Hugo_Symbol=='TP53') %>% select(SubjectID) %>% distinct()

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  filter(SubjectID %in% id_abi) %>%
  mutate(tp53=if_else(SubjectID %in% c1d1.samples$SubjectID, 'mutant','wildtype'))
fit <- survfit(Surv(time, stat) ~ tp53, data = m1)
ggsurvplot(fit, pval = T)

# TP53 野生型和突变型 生存曲线(C1D1期)，只考虑pathogenic, likly pathognic

c1d1.samples <- df %>% 
  filter(StudyVisit=='C1D1') %>%
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, Clinvar,
         CLIN_SIG, time, stat) %>% 
  filter(Hugo_Symbol=='TP53') %>% 
  filter(str_detect(Clinvar, 'Pathogenic|pathogenic')) %>% 
  select(SubjectID) %>% distinct()



# AR mutation----------------------------------------------------------------------

count_sample_gene(long, 'AR')
# 统计long/short中AR突变情况
gene_mutation_count(short, 'AR')

# 统计AR 各mutation的突变情况, 在多少sample中检测到以及在C1期检测到
df %>% filter(Hugo_Symbol=='AR') %>% 
  select(SubjectID, HGVSp_Short) %>% 
  distinct() %>% 
  group_by(HGVSp_Short) %>% count() %>% arrange(desc(n))
  
df %>% filter(Hugo_Symbol=='AR') %>% 
  select(SubjectID, ExternalID, StudyVisit,Variant_Value, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short, StudyVisit) %>% 
  summarise(n = n()) %>% arrange(n) #共24个sample 
df %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, ExternalID, Variant_Value, HGVSp_Short) %>% 
  distinct() %>% group_by(HGVSp_Short) %>% 
  summarise(n = n()) %>% arrange(n) #15个subjectID

# AR 野生型和突变型 生存曲线(C1D1期)
c1d1.samples <- df %>% filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, time, stat) %>% 
  filter(Hugo_Symbol=='AR') %>% select(SubjectID) %>% distinct()

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  filter(SubjectID %in% id_enz) %>%
  mutate(ar=if_else(SubjectID %in% c1d1.samples$SubjectID, 'mutant','wild'))
fit <- survfit(Surv(time, stat) ~ ar, data = m1)
ggsurvplot(fit, pval = T)

# 所有AR mutation + CNV + ARV7
id_ar_cnv <- (df2 %>% filter(Hugo_Symbol == 'AR') %>% 
  filter(StudyVisit == 'C1D1') %>% 
  select(SubjectID) %>% distinct())$SubjectID

id_ar_arv7 <- (df3 %>% filter(StudyVisit == 'C1D1') %>% 
  filter(NameAST == 'AR-V7') %>% 
  select(SubjectID) %>% distinct())$SubjectID

ar_all_id <- unique(c(c1d1.samples$SubjectID, id_ar_cnv, id_ar_arv7))

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  # filter(SubjectID %in% id_enz) %>%
  mutate(AR=if_else(SubjectID %in% ar_all_id, 'mutant','wild'))
fit <- survfit(Surv(time, stat) ~ AR, data = m1)
ggsurvplot(fit, pval = T)



# 对于ARcnv而言-------------------------------------
sort(table(cfdna_cnv$Hugo_Symbol)) #发生CNV最多次的为AR,PIK3CA,PTEN,PPP2R2A
df2 <- cfdna_cnv %>% left_join(d, by = c('SubjectID'='pt. ID')) %>% 
  mutate(StudyVisit=if_else(StudyVisit=='End of Treatment', 'EOT', StudyVisit))
# long/short中CNV出现的频率
df2 %>% filter(SubjectID %in% short) %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(SubjectID) %>% count()
# AR CNV 在不同时期的出现次数
df2 %>% 
  filter(Hugo_Symbol=='AR') %>% 
  select(SubjectID,StudyVisit,Variant_Value) %>% 
  pivot_wider(names_from = StudyVisit, values_from = Variant_Value) %>% 
  mutate(across(where(is.double), ~if_else(is.na(.), 'NA', 'YES')))
# 
hxd2 <- df2 %>% 
  filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID,StudyVisit,Variant_Value) %>% 
  group_by(SubjectID) %>% count()
# plot
df2 %>%
  filter(Hugo_Symbol=='AR') %>% 
  ggplot(., aes(x=StudyVisit, y=time_d)) +
  geom_jitter(width = 0.1) + 
  theme_bw()
  ggrepel::geom_label_repel(aes(label=SubjectID))
# 首次出现AR CNV的样本时间
yy <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  group_by(SubjectID) %>% 
  slice_min(StudyVisit) %>% 
  select(SubjectID,StudyVisit,Variant_Value,time, time_d, stat)
yy$StudyVisit <- factor(yy$StudyVisit, 
                        levels = c('C1D1','C4D1','C7D1','C10D1','C16D1','EOT'))
ggplot(yy, aes(x=StudyVisit, y=time)) +
  geom_jitter(width = 0.1) + 
  theme_bw()

# 分时期对某些基因CNV 进行展示
i_gene <- Hmisc::Cs(AR,PPP2R2A,PTEN,PIK3CA,RB1,MYC,TP53,KDM6A)
i_cnv <- df2 %>% filter(StudyVisit=='C4D1') %>% 
  filter(Hugo_Symbol %in% i_gene) %>% 
  arrange(time) %>% 
  select(SubjectID, Hugo_Symbol, Variant_Value) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Value) %>% 
  select(sort(names(.)))
# i_cnv <- i_cnv[, order(names(i_cnv))]
ff <- psa_meta %>% select(SubjectID,
                          `PSA Percent Change from Baseline`) %>% 
  distinct() %>% 
  arrange(`PSA Percent Change from Baseline`) %>% 
  left_join(i_cnv, by = c('SubjectID'))
M <- ff %>% 
  select(-`PSA Percent Change from Baseline`) %>% 
  column_to_rownames('SubjectID')
ha <- columnAnnotation(bar = anno_barplot(ff$`PSA Percent Change from Baseline`))

# C1D1期AR CNV 出现与否的survival
cnv.s <- function(gene){
  c1d1 <- df2 %>% filter(StudyVisit=='C1D1') %>% 
    filter(Hugo_Symbol==gene) %>% 
    select(SubjectID) %>% distinct()
  c1d1$SubjectID
}
c1d1.ar <- cnv.s('AR')
c1d1.ppp <- cnv.s('PPP2R2A')
c1d1.pten <- cnv.s('PTEN')
c1d1.pik3ca <- cnv.s('PIK3CA')
c1d1.RB1 <- cnv.s('RB1')
c1d1.myc <- cnv.s('MYC')
c1d1.TP53 <- cnv.s('TP53')
c1d1.KDM6A <- cnv.s('KDM6A')
c1d1.CDK6 <- cnv.s('CDK6')
c1d1.CDH1 <- cnv.s('CDH1')
c1d1.atm <- cnv.s('ATM')
c1d1.chek2 <- cnv.s('CHEK2')
c1d1.brca2 <- cnv.s('BRCA2')

all_cnv <- unique(c(c1d1.ar,c1d1.pten,c1d1.pik3ca,c1d1.RB1,c1d1.myc,
                    c1d1.TP53
                    ))
# all_cnv <- unique(c(c1d1.TP53,c1d1.myc,c1d1.atm,c1d1.chek2))

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  filter(SubjectID %in% id_enz) %>%
  mutate(cnv=if_else(SubjectID %in% c1d1.ar, 'mutant','wildtype'))
fit <- survfit(Surv(time, stat) ~ cnv, data = m1)
ggsurvplot(fit, pval = T)


m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  filter(SubjectID %in% id_enz) %>% 
  mutate(ar=if_else(SubjectID %in% c1d1.ar, 2,1)) %>%
  mutate(PPP2R2A=if_else(SubjectID %in% c1d1.ppp, 2,1)) %>% 
  mutate(PTEN=if_else(SubjectID %in% c1d1.pten, 2,1)) %>% 
  mutate(PIK3CA=if_else(SubjectID %in% c1d1.pik3ca, 2,1)) %>% 
  mutate(RB1=if_else(SubjectID %in% c1d1.RB1, 2,1)) %>% 
  mutate(MYC=if_else(SubjectID %in% c1d1.myc, 2,1)) %>% 
  mutate(TP53=if_else(SubjectID %in% c1d1.TP53, 2,1))
  mutate(KDM6A=if_else(SubjectID %in% c1d1.KDM6A, 2,1)) %>% 
  mutate(CDK6=if_else(SubjectID %in% c1d1.CDK6, 2,1)) %>% 
  mutate(CDH1=if_else(SubjectID %in% c1d1.CDH1, 2,1))

reg <- coxph(Surv(time, stat) ~ ar + PPP2R2A+PTEN+PIK3CA+RB1+MYC +TP53,
  #+KDM6A+CDK6+CDH1, 
             data = m1)
summary(reg)
survminer::ggforest(reg, data = m1)


# 对于AR splicing而言----------------------------------------
df3 <- cfrna_splicing %>% left_join(d, by = c('SubjectID'='pt. ID')) %>% 
  mutate(StudyVisit=if_else(StudyVisit=='End of Treatment', 'EOT', StudyVisit))
sort(table(df3$NameAST))
df3 %>% group_by(NameAST) %>% count() %>% arrange(desc(n))
# long/short中出现splicing的频率
df3 %>% filter(SubjectID %in% short) %>% 
  group_by(SubjectID) %>% count()

# AR splicing的分析以AR-V7为主
# 首先还是出现的频率
df3 %>% filter(NameAST == 'AR-V7') %>% filter(SubjectID %in% short) %>% 
  distinct(SubjectID)
df3 %>% filter(NameAST == 'AR-V7') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID) %>% distinct()

# 发生ARV7的样本，其的变异情况
arv7 <- df3 %>% 
  filter(NameAST %in% c('AR-V7'
                        # 'AR-V3', 'AR-V9'
                        )) %>% 
  filter(SupportReads >= 3) %>% 
  select(SubjectID, ExternalID, StudyVisit, NameAST, SupportReads, time, stat) %>%
  arrange(desc(time))

arFL <- df3 %>% filter(NameAST=='AR-FL-7_8') %>% 
  select(SubjectID, StudyVisit, ExternalID, SupportReads)
arFL %>% pivot_wider(values_from = SupportReads, names_from = StudyVisit)

arfl.c1 <- arFL %>% filter(StudyVisit=='C1D1')

aa <- arv7 %>% left_join(arFL, by = c('ExternalID'))

ar_filter <- arv7 %>% filter(time > 0)# >130 or <130
aaid <- ar_filter$ExternalID

(df_arv7 <- df %>% filter(ExternalID %in% aaid) %>% 
  select(SubjectID, ExternalID, StudyVisit, Hugo_Symbol, time) %>% 
  arrange(SubjectID))
# 发生ARV7 同时发生AR mutation的情况
aa <- df_arv7 %>% data.frame() %>% filter(Hugo_Symbol=='AR')
write_delim(aa, 'D:/kintor/arv7.txt', delim = '\t')

# 发生ARV7到事件发生的数据，第一次发生的时间
vv <- df3 %>% filter(NameAST=='AR-V7') %>% 
  group_by(SubjectID) %>% 
  slice_min(StudyVisit) %>% 
  select(SubjectID, StudyVisit, SupportReads, time, time_d) %>% 
  arrange(time_d) %>% 
  filter(SubjectID %in% id_enz)
vv$StudyVisit <- factor(vv$StudyVisit, 
                        levels = c('C1D1','C4D1','C7D1','C16D1','EOT'))
ggplot(vv, aes(x=StudyVisit, y=time_d)) +
  geom_jitter(width = 0.1) + 
  theme_bw() +
  ggrepel::geom_label_repel(aes(label=SubjectID))


# plot --------------------------------------------------------------------

# bar plot
aa <- df %>% select(SubjectID, drug, time) %>% distinct() %>% arrange(drug,time)
  filter(SubjectID %in% id_enz)
aa$SubjectID <- factor(aa$SubjectID, levels = aa$SubjectID)

ggplot(aa) + geom_bar(aes(x = SubjectID, y = time,
                          fill = drug
                          ), stat = 'identity') + theme_bw() +
  ggsci::scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = c(0.1,0.8)
        )


# heatmap 
library(ComplexHeatmap)

# C1D1期某个基因的突变频率
filter_gene <- function(gene, visit='C1D1'){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}
mut_tp53 <- filter_gene('TP53')
mut_ar <- filter_gene('AR')
mut_atm <- filter_gene('ATM')
mut_ctnnb1 <- filter_gene('CTNNB1')

ar_cnv <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')

final_data <- psa_meta %>% select(SubjectID,`PSA Percent Change from Baseline`) %>% 
  distinct() %>% 
  arrange(`PSA Percent Change from Baseline`) %>% 
  left_join(mut_tp53, by = c('SubjectID')) %>% 
  left_join(mut_atm, by = c('SubjectID')) %>%
  left_join(mut_ctnnb1, by = c('SubjectID')) %>%
  left_join(mut_ar, by = c('SubjectID')) %>%
  left_join(ar_cnv, by = c('SubjectID'))

M <- final_data %>% 
  select(-`PSA Percent Change from Baseline`) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$`PSA Percent Change from Baseline`))

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
        width = ncol(M)*unit(60, "mm"),
        height = nrow(M)*unit(1, "mm"),
        # top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)

# 选择某个基因的mutation频率
filter_aa <- function(gene, aa){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit=='C1D1') %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}
ar_L702H <- filter_aa('AR','p.L702H')

ar_gene <- c('p.H875Y','p.T878A','p.L702H','p.T878S',
             'p.R787Q'
)
final_data <- psa_meta %>% select(SubjectID,
                                  `PSA Percent Change from Baseline`) %>% 
  distinct() %>% 
  arrange(`PSA Percent Change from Baseline`)
for(i in ar_gene){
  foo <- filter_aa('AR', i)
  final_data <- final_data %>%
    left_join(foo, by = c('SubjectID'))
  
}

ar_v7 <- df3 %>% filter(NameAST=='AR-V7') %>% 
  filter(StudyVisit=='EOT') %>% 
  select(SubjectID, SupportReads) %>% 
  rename(arv7=SupportReads)

ar_cnv <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='EOT') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')

final_data %<>% left_join(ar_cnv, by = c('SubjectID')) %>% 
  left_join(ar_v7, 'SubjectID')

M <- final_data %>% 
  select(-`PSA Percent Change from Baseline`) %>%
  column_to_rownames('SubjectID')

ha <- columnAnnotation(bar = anno_barplot(final_data$`PSA Percent Change from Baseline`))



#
# re-analysis of cfdna_SNV ------------------------------------------------

cal_gene_appear <- function(gene){
  foo <- df %>% filter(Hugo_Symbol==gene) %>% group_by(ExternalID) %>% 
    count()
  return(dim(foo)[1])
}
cal_gene_appear('TP53')

y <- vector()
for (i in unique(df$Hugo_Symbol)) {
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
m1 <- df %>% select(SubjectID, time) %>% distinct()
cfrna_exp <- readxl::read_excel('D:/kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
                                sheet = 'RNA_Expression') %>% 
  left_join(m1, by='SubjectID') %>% 
  mutate(StudyVisit=if_else(StudyVisit=='End of Treatment', 'EOT', StudyVisit))

library(ComplexHeatmap)

fpkm <- cfrna_exp %>% select(-c(SubjectID,StudyVisit,Gene,time,time_d)) %>% 
  column_to_rownames('ExternalID') %>% t() %>% 
  as.data.frame()
cc <- cfrna_exp %>% mutate(class = if_else(SubjectID %in% short, 'short', 'long'))
anno <- cc %>% select(ExternalID, class) %>% 
  column_to_rownames('ExternalID')

p_heat <- pheatmap::pheatmap(
  log10(fpkm + 1),
  cluster_cols = T,
  cluster_rows = T,
  show_rownames = F,
  show_colnames = F,
  annotation_col = anno
)
# use complexheatmap
col_anno <- cc %>% select(ExternalID, time_d) %>% 
  column_to_rownames('ExternalID')
col_anno %<>% replace(is.na(.),0)
ha <- columnAnnotation(time = anno_barplot(col_anno))
Heatmap(
  log10(fpkm + 1),
  name = 'log10(fpkm)',
  top_annotation = ha,
  column_names_gp = gpar(fontsize=8),
  row_names_gp = gpar(fontsize=1),
  show_column_names = F
)


# only for C1D1期
fpkm_c1d1 <- cfrna_exp %>% filter(StudyVisit=='C1D1') %>% 
  select(-c(ExternalID,StudyVisit,Gene,time)) %>% 
  column_to_rownames('SubjectID') %>% t() %>% 
  as.data.frame()

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
  cluster_columns = T
)

# 只选择感兴趣的基因，
i_gene <- Hmisc::Cs(TP53,AR,CHEK2,ATM,PIK3CA)
s <- cfrna_exp %>% 
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

cfrna_snv_others <- readxl::read_excel('D:/us1002/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_20210508.xlsx',
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
              unite('key', c(ExternalID,Hugo_Symbol,HGVSp), remove = F))$key %>% 
  unique()

rna_other_key <- (cfrna_snv_others %>% 
                    unite('key', c(ExternalID,GeneSymbol,HGVSp), remove = F))$key %>% 
  unique()

sum(rna_other_key %in% dna_key)
sum(rna_key %in% dna_key)
print(length(rna_key))
print(length(dna_key))
sum(rna_key %in% rna_other_key)




# PSA ---------------------------------------------------------------------

d_id <- unique(d$`pt. ID`)
id_44 <- unique(df$SubjectID)

psa_f <- readxl::read_excel('D:/us1002/l_psa_reduction from PHoebe 20201215.xlsx',
                            sheet = 'l_psa_reduction'
                            )

psa_max <- psa_f %>% 
  filter(`Patient Number` %in% id_44) %>% 
  # filter(!grepl('Unscheduled', Visit)) %>% 
  filter(!str_detect(Visit,'Unscheduled')) %>% 
  filter(Visit != 'Screening') %>%
  group_by(`Patient Number`) %>% 
  slice_max(abs(`PSA Percent Change from Baseline`)) %>% 
  arrange(`PSA Percent Change from Baseline`)
# write_delim(psa_max, 'D:/kintor/psa_max.txt', delim = '\t')
psa_max <- read_delim('D:/kintor/psa_max.txt', 
                      delim = '\t')

psa_max$`Patient Number`<- factor(psa_max$`Patient Number`, 
                                  levels = psa_max$`Patient Number`)


ggplot(psa_max) + geom_bar(aes(x = `Patient Number`, 
                               y = `PSA Percent Change from Baseline`/100), 
                      stat = 'identity')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = 'none'
  ) + ggsci::scale_fill_aaas() +
  scale_y_break(c(1,5), scales = 0.2, ticklabels = NULL) +
  ylab('Change from Baseline')

#
psa_meta <- d %>% filter(`pt. ID` %in% id_44) %>% 
  left_join(psa_max, by = c('pt. ID'='Patient Number')) %>% 
  arrange(`PSA Percent Change from Baseline`) %>% 
  rename('SubjectID'=`pt. ID`) %>% 
  mutate(`PSA Percent Change from Baseline`=`PSA Percent Change from Baseline`/100) %>% 
  arrange(`PSA Percent Change from Baseline`) %>% 
  mutate(class = if_else(`PSA Percent Change from Baseline`>0,'up','down')) %>% 
  mutate(`PSA Percent Change from Baseline`=if_else(`PSA Percent Change from Baseline`>1,1,
                 `PSA Percent Change from Baseline`))
psa_meta$SubjectID <- factor(psa_meta$SubjectID, levels = psa_meta$SubjectID)
# write_delim(psa_meta, 'D:/kintor/psa_meta.txt', delim = '\t')

ggplot(psa_meta) + geom_bar(aes(x = `SubjectID`, 
                               y = `PSA Percent Change from Baseline`,
                               fill = class
                               ), 
                           stat = 'identity')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = 'none'
  ) + ggsci::scale_fill_jama() +
  scale_y_break(c(1,5), scales = 0.2, ticklabels = NULL) +
  ylab('Change from Baseline')
  


#
# cfdna concentration -----------------------------------------------------
# 浓度变化的比例, 或者看一下浓度和time之间的相关性
# sperman correlation

haha <- cfdna_qc %>% filter(StudyVisit=='C1D1') %>% 
  arrange(SubjectID)

hehe <- t %>% arrange(SubjectID)

heha <- bind_cols(haha$DNAYield, hehe$time)
names(heha) <- c('DNA','time')
correlation::cor_test(heha, x = 'DNA', y = 'time', 
                      method = 'spearman')

# 定义mutation burden


# 浓度和生存的关系


# 关于mutation 突变频率






# AR mutation & CNV  & TP53--------------------------------------------------

c1d1.samples <- df %>% 
  filter(StudyVisit=='C1D1') %>%
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, time, stat) %>% 
  filter(Hugo_Symbol=='AR') %>% select(SubjectID) %>% 
  distinct() %>% pull(SubjectID)
c1d1.tp53 <- df %>% 
  filter(StudyVisit=='C1D1') %>%
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, time, stat) %>% 
  filter(Hugo_Symbol=='TP53') %>% select(SubjectID) %>% 
  distinct()
cnv.s <- function(gene){
  c1d1 <- df2 %>% 
    filter(StudyVisit=='C1D1') %>%
    filter(Hugo_Symbol==gene) %>% 
    select(SubjectID) %>% distinct()
  c1d1$SubjectID
}
c1d1.ar <- cnv.s('AR')

v7.s <- df3 %>% filter(NameAST=='AR-V7') %>% 
  filter(StudyVisit=='C1D1') %>% 
  select(SubjectID) %>% distinct()

c1d1.samples <- base::union(c1d1.ar, c1d1.samples)
c1d1.samples <- base::union(c1d1.samples, v7.s$SubjectID)
# c1d1.samples <- base::union(c1d1.samples, c1d1.tp53$SubjectID)
length(unique(c1d1.samples))

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  filter(SubjectID %in% id_enz) %>% 
  mutate(ar=if_else(SubjectID %in% c1d1.samples, 'mutant','wildtype'))

fit <- survfit(Surv(time, stat) ~ ar, data = m1)
ggsurvplot(fit, pval = T)

# 
ggsurv <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = m1,               # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Male", "Female")    # change legend labels.
)
ggsurv



# 重新定义的EOT ----------------------------------------------------------------
# 对EOT期没有检测的样本进行定义
eot.na <- cfdna_qc %>% select(SubjectID, StudyVisit, LibraryYield) %>% 
  pivot_wider(names_from = StudyVisit, values_from = LibraryYield)

glimpse(eot.na)

foo <- eot.na %>% 
  mutate(C4D1 = ifelse(is.na(C4D1),C4D1, 'C4D1')) %>% 
  mutate(C7D1 = ifelse(is.na(C7D1),C7D1, 'C7D1'),
         C10D1 = ifelse(is.na(C10D1),C10D1, 'C10D1'),
         C13D1 = ifelse(is.na(C13D1),C13D1, 'C13D1'),
         C16D1 = ifelse(is.na(C16D1),C16D1, 'C16D1'),
         C19D1 = ifelse(is.na(C19D1),C19D1, 'C19D1')
         ) %>% 
  mutate(EOTT = ifelse(is.na(EOT),
                       ifelse(is.na(C19D1),ifelse(
                         is.na(C16D1),ifelse(
                           is.na(C13D1),ifelse(
                             is.na(C10D1),ifelse(
                               is.na(C7D1),ifelse(
                                 is.na(C4D1),NA,C4D1
                               ),C7D1
                             ),C10D1
                           ),C13D1
                         ),C16D1
                       ),C19D1)
                       ,EOT)) %>% 
  filter(str_detect(EOTT, 'C'))

id_c19 <- (foo %>% filter(EOTT=='C19D1'))$SubjectID
id_c16 <- (foo %>% filter(EOTT=='C16D1'))$SubjectID
id_c13 <- (foo %>% filter(EOTT=='C13D1'))$SubjectID
id_c10 <- (foo %>% filter(EOTT=='C10D1'))$SubjectID
id_c7 <- (foo %>% filter(EOTT=='C7D1'))$SubjectID
id_c4 <- (foo %>% filter(EOTT=='C4D1'))$SubjectID
  
filter_gene <- function(gene, visit='EOT'){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    ungroup() %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}

filter_gene_id <- function(gene, id, visit='C19D1'){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    ungroup() %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}

aa <- function(gene){
  eot22 <- filter_gene(gene)
  f19 <- filter_gene_id(gene, id_c19, visit = 'C19D1')
  f16 <- filter_gene_id(gene, id_c16, visit = 'C16D1')
  f13 <- filter_gene_id(gene, id_c13, visit = 'C13D1')
  f10 <- filter_gene_id(gene, id_c10, visit = 'C10D1')
  f7 <- filter_gene_id(gene, id_c7, visit = 'C7D1')
  f4 <- filter_gene_id(gene, id_c4, visit = 'C4D1')
  ff <- eot22 %>% 
    bind_rows(f19) %>% 
    bind_rows(f16) %>% 
    bind_rows(f13) %>% 
    bind_rows(f10) %>% 
    bind_rows(f7) %>% 
    bind_rows(f4)
  ff
}


mut_tp53 <- aa('TP53')
mut_ar <- aa('AR')
mut_atm <- aa('ATM')
mut_ctnnb1 <- aa('CTNNB1')

ar_cnv_e <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='EOT') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')

filter.arcnv <- function(gene,visit, id){
  df2 %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    select(SubjectID,Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    rename('AR_CNV'='Variant_Value')
}

bb <- function(gene){
  f19 <- filter.arcnv(gene, id_c19, visit = 'C19D1')
  f16 <- filter.arcnv(gene, id_c16, visit = 'C16D1')
  f13 <- filter.arcnv(gene, id_c13, visit = 'C13D1')
  f10 <- filter.arcnv(gene, id_c10, visit = 'C10D1')
  f7 <- filter.arcnv(gene, id_c7, visit = 'C7D1')
  f4 <- filter.arcnv(gene, id_c4, visit = 'C4D1')
  ff <- ar_cnv_e %>% 
    bind_rows(f19) %>% 
    bind_rows(f16) %>% 
    bind_rows(f13) %>% 
    bind_rows(f10) %>% 
    bind_rows(f7) %>% 
    bind_rows(f4)
  ff
}

ar_cnv <- bb('AR')

final_data <- df %>% select(SubjectID, time) %>% 
  distinct() %>% 
  arrange(time) %>% 
  # filter(SubjectID %in% id_enz) %>% 
  left_join(mut_tp53, by = c('SubjectID')) %>% 
  left_join(mut_atm, by = c('SubjectID')) %>%
  left_join(mut_ctnnb1, by = c('SubjectID')) %>%
  left_join(mut_ar, by = c('SubjectID')) %>%
  left_join(ar_cnv, by = c('SubjectID'))

M <- final_data %>% 
  select(-time) %>%
  column_to_rownames('SubjectID')

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
        col = circlize::colorRamp2(c(0,5), c("white", "firebrick3")),
        width = ncol(M)*unit(40, "mm"),
        height = nrow(M)*unit(2, "mm"),
        # top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)

# 选择某个基因的mutation频率
filter_aa <- function(gene, aa){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit=='EOT') %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}

filter_aa_id <- function(gene, aa, visit, id){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>%
    filter(HGVSp_Short==aa) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    filter(SubjectID %in% id) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    rename((!!aa) := Variant_Value)
  res
}

ar_gene <- c('p.H875Y','p.T878A','p.L702H','p.T878S',
             'p.R787Q'
)
final_data_eot <- df %>% select(SubjectID,
                                  time) %>% 
  distinct() %>% 
  arrange(time) %>% 
  filter(!SubjectID %in% c(id_c19,id_c16,id_c13,id_c10,id_c7,id_c4))
for(i in ar_gene){
  foo <- filter_aa('AR', i)
  final_data_eot <- final_data_eot %>%
    left_join(foo, by = c('SubjectID'))
}

f.loop <- function(id, visit){
  f_data <- df %>% select(SubjectID,
                          time) %>% 
    distinct() %>% 
    arrange(time) %>% 
    filter(SubjectID %in% id)
  for(i in ar_gene){
    foo <- filter_aa_id('AR', i, visit = visit, id = id)
    f_data <- f_data %>%
      left_join(foo, by = c('SubjectID'))
  }
  f_data
  }

f19 <- f.loop(id_c19, visit = 'C19D1')
f16 <- f.loop(id_c16, visit = 'C16D1')
f13 <- f.loop(id_c13, visit = 'C13D1')
f10 <- f.loop(id_c10, visit = 'C10D1')
f7 <- f.loop(id_c7, visit = 'C7D1')
f4 <- f.loop(id_c4, visit = 'C4D1')

final_data <- final_data_eot %>% 
  bind_rows(f19) %>% 
  bind_rows(f16) %>% 
  bind_rows(f13) %>% 
  bind_rows(f10) %>% 
  bind_rows(f7) %>% 
  bind_rows(f4) %>% 
  arrange(time)


# lrm for psa -------------------------------------------------------------
# 一个较为完备的线性回归
library(naniar)

visdat::vis_dat(psa_meta)
psa_meta %>% miss_var_summary() # statistic NA value











# 其他mutation的分布 -----------------------------------------------------------
# 全部mutation再C1D1期的分布
c1d1.all.mutation <- df %>% filter(StudyVisit == 'C1D1') %>% 
  filter(Variant_Value >= 0) %>% 
  group_by(Hugo_Symbol) %>% 
  count() %>% arrange(desc(n))

df %>% filter(StudyVisit == 'C1D1') %>% 
  select(SubjectID, Hugo_Symbol) %>% distinct() %>% 
  group_by(Hugo_Symbol) %>% count() %>% arrange(desc(n))

# write_delim(c1d1.all.mutation, 'D:/us1002/c1d1.all.mutation.txt', delim = '\t')

c1d1.all.cnv <- df2 %>% 
  filter(StudyVisit=='C1D1') %>% 
  group_by(Hugo_Symbol) %>% 
  count() %>% arrange(desc(n))

# write_delim(c1d1.all.cnv, 'D:/us1002/c1d1.all.cnv.txt', delim = '\t')

# 某基因在C1D1期出现与否的生存分析
# TP53,DNMT3A,TET2,ATM,ASXL1,CHEK2,AR,CBL
c1d1.samples <- df %>% 
  filter(StudyVisit=='C1D1') %>%
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, Clinvar,
         CLIN_SIG, time, stat) %>% 
  filter(Hugo_Symbol=='CHEK2') %>% 
  filter(str_detect(Clinvar, 'Pathogenic|pathogenic')) %>%
  select(SubjectID) %>% distinct() %>% 
  pull(SubjectID)
print(length(c1d1.samples))

m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  # filter(SubjectID %in% id_abi) %>%
  mutate(CHEK2=if_else(SubjectID %in% c1d1.samples, 'mutant','wildtype'))

fit <- survfit(Surv(time, stat) ~ CHEK2, data = m1)
ggsurvplot(fit, pval = T)



# TP53 通路相关基因, C1D1期生存分析, mutation+CNV
tp53_gene <- Hmisc::Cs(TP53, ATM, CHEK2, MYC, RB1)

c1d1.samples <- df %>% filter(StudyVisit=='C1D1') %>% 
  select(SubjectID, Variant_Value, Clinvar, Hugo_Symbol,HGVSp_Short, time, stat) %>% 
  filter(Hugo_Symbol %in% tp53_gene) %>% 
  # filter(str_detect(Clinvar, 'Pathogenic|pathogenic')) %>% 
  select(SubjectID) %>% distinct() %>% 
  pull(SubjectID)

cnv.s <- function(gene){
  c1d1 <- df2 %>% filter(StudyVisit=='C1D1') %>% 
    filter(Hugo_Symbol %in% gene) %>% 
    select(SubjectID) %>% distinct()
  c1d1$SubjectID
}

cnv_id <- cnv.s(tp53_gene)

all_sample <- unique(c1d1.samples, cnv_id)
print(length(all_sample))


m1 <- df %>% select(SubjectID, time, stat) %>% distinct() %>%
  mutate(TP53=if_else(SubjectID %in% all_sample, 'mutant','wildtype'))
fit <- survfit(Surv(time, stat) ~ TP53, data = m1)
ggsurvplot(fit, pval = T)


#
# delet 6 patients --------------------------------------------------------

delet6 <- c(
  '006-019',
  '002-005',
  '006-018',
  '002-004',
  '005-012',
  '006-017'
)

cfdna_qc %>% filter(!SubjectID %in% delet6) %>% 
  select(ExternalID) %>% distinct()

cfdna_qc %>% filter(!SubjectID %in% delet6) %>% 
  select(SubjectID) %>% distinct()

psa_max %>% 
  filter(`PSA Percent Change from Baseline` < -50) %>% 
  filter(!`Patient Number` %in% delet6)

tt <- t %>% filter(!SubjectID %in% delet6)

df %>% filter(Hugo_Symbol=='AR') %>% 
  # filter(StudyVisit=='C1D1') %>%
  filter(!SubjectID %in% delet6) %>% 
  filter(HGVSp_Short=='p.T878S') %>%
  select(SubjectID) %>% distinct()

hxd2 <- df2 %>% 
  filter(Hugo_Symbol=='AR') %>% 
  # filter(StudyVisit=='C1D1') %>% 
  filter(!SubjectID %in% delet6) %>% 
  select(SubjectID,StudyVisit,Variant_Value) %>% 
  group_by(SubjectID) %>% count()

length(base::union(hxd$SubjectID, hxd2$SubjectID))

df3 %>% filter(NameAST == 'AR-V7') %>% 
  filter(StudyVisit == 'C1D1') %>%
  filter(!SubjectID %in% delet6) %>% 
  select(SubjectID) %>% distinct()

df %>% 
  filter(StudyVisit=='C1D1') %>%
  select(SubjectID, Variant_Value, Hugo_Symbol,HGVSp_Short, Clinvar,
         CLIN_SIG, time, stat) %>% 
  filter(!SubjectID %in% delet6) %>% 
  filter(Hugo_Symbol=='TP53') %>% 
  filter(str_detect(Clinvar, 'Pathogenic|pathogenic')) %>% 
  select(SubjectID) %>% distinct()







# 按照不同用药前分组 ---------------------------------------------------------------

id_abi <- (df %>% filter(drug=='abiraterone') %>% select(SubjectID) %>% 
  distinct())$SubjectID

id_enz <- (df %>% filter(drug=='Enzalutamide') %>% select(SubjectID) %>% 
             distinct())$SubjectID


filter_gene <- function(gene, visit='C4D1'){
  res <- df %>% filter(Hugo_Symbol==gene) %>% 
    filter(StudyVisit==visit) %>% 
    arrange(SubjectID, Variant_Value) %>% 
    select(SubjectID,Variant_Value) %>% 
    group_by(SubjectID) %>% 
    slice_max(Variant_Value) %>% 
    dplyr::rename(!!sym(gene) := Variant_Value)
  res
}
mut_tp53 <- filter_gene('TP53')
mut_ar <- filter_gene('AR')
mut_atm <- filter_gene('ATM')
mut_ctnnb1 <- filter_gene('CTNNB1')

ar_cnv <- df2 %>% filter(Hugo_Symbol=='AR') %>% 
  filter(StudyVisit=='C4D1') %>% 
  select(SubjectID, Variant_Value) %>% 
  rename('AR_CNV'='Variant_Value')

final_data <- df %>% select(SubjectID, time) %>% 
  distinct() %>% 
  arrange(time) %>% 
  # filter(SubjectID %in% id_enz) %>% 
  left_join(mut_tp53, by = c('SubjectID')) %>% 
  left_join(mut_atm, by = c('SubjectID')) %>%
  left_join(mut_ctnnb1, by = c('SubjectID')) %>%
  left_join(mut_ar, by = c('SubjectID')) %>%
  left_join(ar_cnv, by = c('SubjectID'))

M <- final_data %>% 
  select(-time) %>%
  column_to_rownames('SubjectID')

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
        col = circlize::colorRamp2(c(0,5), c("white", "firebrick3")),
        width = ncol(M)*unit(40, "mm"),
        height = nrow(M)*unit(2, "mm"),
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = F,
        # column_names_gp = gpar(fontsize=1)
)









# 统计同时发生CNV 和 mutation的样本
find_overlap <- function(gene, study='C1D1', ...){
  cnv_id <- df2 %>% 
    filter(StudyVisit==study) %>% 
    filter(Hugo_Symbol==gene) %>% 
    pull(SubjectID)
  
  df %>% filter(SubjectID %in% cnv_id) %>% 
    filter(Hugo_Symbol == gene) %>% 
    pull(SubjectID) %>% unique()
}

find_overlap('AR')

