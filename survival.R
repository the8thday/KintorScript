

library(survival)
library(survminer)
library(tidyverse)
library(lubridate)
library(magrittr)
# library(forestplot)

cfdna_snv <- readxl::read_excel('/Users/congliu/OneDrive/kintor/GT0918-US-1002_cumulative_cfDNA_cfRNA_data_202010120.xlsx',
                                sheet = 'cfDNA_SNV')
metadata <- readxl::read_excel('D:/kintor/GT0918-US-1002_Subject_status_2020-1014-general_dates-TH_12Mar2021.xlsx')
d <- metadata %>% 
  mutate(time = `EOT`-`Date of 1st Dosing`) %>% 
  mutate(stat = if_else(`Treatment discontunation`%in% 
                          c('Bone or Radiographic Disease Progression')|
                          `if others, specify`%in% c('Rise in PSA',
                                                     'PI decision due to rising PSA',
                                                     'rising PSA and intolerable toxicities'
                                                     ), 
                        2, 1)) %>% 
  select(c('pt. ID', 'dosing', 'Prior Treat', 'time', 'stat')) %>% 
  rename(drug = `Prior Treat`)
d$time %<>% as.vector

df <- cfdna_snv %>% left_join(d, by = c('SubjectID'='pt. ID'))

df_for_cfdna <- df %>% select(c(SubjectID, ExternalID, StudyVisit, Hugo_Symbol,
                                   Variant_Value, HGVSc, HGVSp, dosing, drug, time, stat
                                   )) %>% 
  distinct() %>% 
  unite('variant', c(Hugo_Symbol, HGVSp), sep = ':', remove = F) %>% 
  # group_by(variant) %>% filter(n()>2) %>% 
  mutate(dosing = ifelse(dosing=='400mg', 1, 2)) %>% 
  mutate(drug = ifelse(drug=='abiraterone', 1, 2))

df_for_all <- d %>% 
  mutate(dosing = ifelse(dosing=='400mg', 1, 2)) %>% 
  mutate(drug = ifelse(drug=='abiraterone', 1, 2))
# 1为400mg
# 1为abiraterone
df_for_all$dosing <- factor(df_for_all$dosing, levels = c(1,2), 
                            labels = c('400mg', '500mg'))
df_for_all$drug <- factor(df_for_all$drug, levels = c(1,2), 
                            labels = c('abiraterone', 'Enzalutamide'))


# cox reg -----------------------------------------------------------------

reg <- coxph(Surv(time, stat) ~ drug, 
             data = df_for_all)
summary(reg)
survminer::ggforest(reg, data = df_for_all)

reg2 <- coxph(Surv(time, stat) ~ drug + dosing, 
              data = df_for_all)
survminer::ggforest(reg2, data = df_for_all)
summary(reg2)
ggcoxadjustedcurves(reg2)
ggcoxdiagnostics(reg2)

# pec, rms, 以及PH等比性 res_cox <- cox.zph(reg2)K
reg3 <- coxph(Surv(time, stat) ~ drug + dosing, 
              data = df_for_all)
z <- cox.zph(reg2)
ggcoxzph(z)

# nomograph
library(rms)

df_for_rms <- df_for_all
df_for_rms$dosing <- factor(df_for_rms$dosing,
                            levels = c(1,2),
                            labels = c('400mg', '500mg')
                            )
df_for_rms$drug <- factor(df_for_rms$drug,
                          levels = c(1, 2),
                          labels = c('abiraterone', 'Enzalutamide')
                          )
dd <- datadist(df_for_rms)
options(datadist='dd')

coxm <- cph(Surv(time, stat) ~ drug+dosing,
            x = T,
            y = T,
            surv = T,
            data = df_for_rms
            )
surv <- Survival(coxm)
surv1 <- function(x){surv(100, lp=x)}
surv2 <- function(x)surv(200, lp=x)
med <- Quantile(coxm)

nom <- nomogram(coxm,
                fun = function(x){med(lp=x)}
                )
plot(nom, xfrac=0.3)

# c-index, 补数
summary(reg2)$concordance #即是C-index
S <- Surv(df_for_all$time, df_for_all$stat)
Hmisc::rcorrcens(S~predict(coxm), outx=TRUE)

# 校准曲线
cal <- calibrate(coxm,
                 cmethod = 'KM',
                 method = 'boot',
                 u = 100
                 )


# km plot ---------------------------------------------------------------------
# 1为截断数据，2为事件发生
fit <- survfit(Surv(time, stat) ~ drug, data = d)
summary(fit)
surv_summary(fit)
g1 <- ggsurvplot(fit, pval = T)
# log rank
survdiff1 <- survdiff(Surv(time, stat) ~ drug, data = d)
survdiff1$chisq
pValue <- 1 - pchisq(survdiff1$chisq, length(survdiff1$n) - 1)
surv_pvalue(fit = fit, data = d)
hr <- survdiff1$obs[1] * survdiff1$exp[2]/(survdiff1$obs[2] * survdiff1$exp[1])
plot1 <- ggsurvplot(fit, pval = T)
plot1$plot + annotate()


fit2 <- survfit(Surv(time, stat) ~ dosing, data = d)
g2 <- ggsurvplot(fit2, pval = T)

arrange_ggsurvplots(list(g1, g2))

ggsurvevents(fit = fit, data = d)
surv_cutpoint() # 为连续变量分配最佳cutpoint



# Customized survival curves
ggsurvplot(fit, 
           data = m1,
           surv.median.line = "hv",
           # Change legends: title & labels
           legend.title = "class",
           # legend.labs = c("Male", "Female"),
           censor.shape="|", censor.size = 4,
           # Add p-value and tervals
           pval = TRUE,
           pval.method = TRUE,
           # pval.coord = c(0, 0.03), #调节Pval的位置
           # conf.int = TRUE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           # palette = c("#E7B800", "#2E9FDF"), 
           # censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小
           ylab="Cumulative survival (percentage)",xlab = " Time (Days)", #更改横纵坐标
           palette = 'aaas',
           ggtheme = theme_survminer() +
             theme(legend.text = element_text(size = 15),
                   legend.title = element_text(size = 15)
                   )
)
###添加HR ,CI ,P
res_cox<-coxph(Surv(time, status) ~sex, data=lung)
p3$plot = p3$plot + ggplot2::annotate("text",x = 50, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
  ggplot2::annotate("text",x = 50, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 50, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
# https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
data.survdiff <- survdiff(Surv(time, status) ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))


fit <- survfit(Surv(time, stat) ~ dosing, data = df_for_survival)
ggsurvplot(fit = fit, pval = T)
ggsurvplot_facet(fit, 
                 data = df_for_survival, 
                 facet.by = "StudyVisit",
                 palette = "jco",
                 pval = TRUE)


# survival for tp53 -------------------------------------------------------
# "p.V216M" "p.R175H" "p.R273H"
aa <- as.character(quote(
  c(p.G244D,
    p.G244C,
    p.Y163C,
    p.H179Q,
    p.Y220C,
    p.R196Q,
    p.T125M,
    p.P278L,
    p.F270S,
    p.Y234C,
    p.R175G,
    p.R248Q,
    p.N239D,
    p.E285K,
    p.C275F,
    p.C275Y,
    p.R290H,
    p.R280T,
    p.I254N,
    p.D281_R283delinsG,
    p.C277_G279delinsW
  )
))[-1]
aa <- c(aa, 'p.D352Gfs*29','p.R196*')

aa_130 <- as.character(quote(
  c(p.G244D,
    p.G244C,
    p.Y163C,
    p.H179Q,
    p.R196Q,
    p.T125M,
    p.P278L,
    p.F270S,
    p.Y234C,
    p.R175G,
    p.R248Q,
    p.N239D,
    p.E285K,
    p.C275F
  )
))[-1]
aa_130 <- c(aa_130, 'p.D352Gfs*29','p.R196*')

dnmt3a <- c(
  "p.V895M",
  "p.R771Q",
  "p.E756*",
  "p.L703Wfs*2",
  "p.V897D",
  "p.R887Efs*34",
  "p.W860Cfs*21",
  "p.P700A",
  "p.F827Lfs*4",
  "p.E733*",
  "p.R729G",
  "p.M880V",
  "p.R736C"
)

# 按照SubjectID
df_p53 <- df %>% mutate(stat_p53 = ifelse(Hugo_Symbol=='TP53' & 
                                            HGVSp_Short %in% aa,
                                          'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)


d_tp53 <- d %>% right_join(df_p53, by = c('pt. ID'='SubjectID'))

# 去掉后面三个较为异常的样本
d_tp53 %<>% filter(! `pt. ID` %in% c('008-005','006-003','005-001'))

fit <- survfit(Surv(time, stat) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

# 按照ExternalID
foo <- df %>% select(ExternalID, time, time_d, stat) %>% distinct()
df_p53 <- df %>% mutate(stat_p53 = ifelse(Hugo_Symbol=='TP53' & 
                                            HGVSp_Short %in% aa,
                                          'yes', 'no')) %>%
  # select(SubjectID, stat_p53, ExternalID, time, time_d, stat) %>%
  select(ExternalID, stat_p53) %>%
  distinct() %>% 
  arrange(ExternalID, desc(stat_p53)) %>% 
  group_by(ExternalID) %>% 
  top_n(1) %>% ungroup()
d_tp53 <- df_p53 %>% left_join(foo, by = c('ExternalID'))
fit <- survfit(Surv(time_d, stat) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)

write_delim(d_tp53, 'D:/kintor/d_tp53.txt', delim = '\t')


# 只考虑患者入组时的状态
ruzu <- df %>% filter(StudyVisit=='C1D1')
df_p53 <- ruzu %>% mutate(stat_p53 = ifelse(Hugo_Symbol=='TP53' & 
                                            HGVSp_Short %in% aa,
                                          'yes', 'no')) %>% 
  select(SubjectID, stat_p53) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)

d_tp53 <- d %>% right_join(df_p53, by = c('pt. ID'='SubjectID'))
fit <- survfit(Surv(time, stat) ~ stat_p53, data = d_tp53)
ggsurvplot(fit, pval = T)


# 从ARV7 出现与否的角度进行分析
# 首先是按照ExternalID进行
foo <- df %>% select(ExternalID, time, time_d, stat) %>% distinct()
df_arv7 <- df3 %>% mutate(stat_arv7 = if_else(NameAST=='AR-V7', 'positive', 'negtive')) %>% 
  select(ExternalID, stat_arv7) %>%
  distinct() %>% 
  arrange(ExternalID, desc(stat_arv7)) %>% 
  group_by(ExternalID) %>% 
  top_n(1) %>% ungroup()
d_arv7 <- df_arv7 %>% right_join(foo, by = c('ExternalID'))
d_arv7 %<>% replace(is.na(.),'negtive')
d_arv7 %<>% filter(time_d > 0)
fit <- surv_fit(Surv(time_d, stat) ~ stat_arv7, data = d_arv7)
ggsurvplot(fit, pval = T)

# 其次是按照患者的入组状态，即C1D1期
ruzu <- df3 %>% filter(StudyVisit=='C1D1')
df_arv7 <- ruzu %>% mutate(stat_arv7 = if_else(NameAST=='AR-V7', 'positive', 'negtive')) %>% 
  select(SubjectID, stat_arv7) %>% distinct() %>% 
  arrange(SubjectID) %>% 
  group_by(SubjectID) %>% 
  top_n(1)

d_arv7 <- d %>% right_join(df_arv7, by = c('pt. ID'='SubjectID'))
fit <- survfit(Surv(time, stat) ~ stat_arv7, data = d_arv7)
ggsurvplot(fit, pval = T)




# 批量 ----------------------------------------------------------------------

library(ezcox)

cc <- c(
  "c.731G>A",
  "c.730G>T",
  "c.488A>G",
  "c.537T>A",
  "c.659A>G",
  "c.587G>A",
  "c.374C>T",
  "c.833C>T",
  "c.809T>C",
  "c.701A>G",
  "c.1055_1056del",
  "c.523C>G",
  "c.743G>A",
  "c.715A>G",
  "c.853G>A",
  "c.824G>A",
  "c.869G>A",
  "c.839G>C",
  "c.761T>A",
  "c.842_847del",
  "c.831_836del_del"
)





