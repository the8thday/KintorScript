
# plot --------------------------------------------------------------------


library(ggpubr)


foo <- d %>% select(c(`pt. ID`,dosing, drug, time)) %>% 
  distinct()
p <- ggboxplot(foo, 
          x = 'dosing',
          y = 'time',
          color = 'dosing',
          palette = c("#00AFBB", "#E7B800", "#FC4E07")
)
my_comparisons <- list( c("400mg", "500mg"))
p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) 

# pathview ----------------------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(KEGGREST)

genelist = as.character(quote(c(ABL1,AKT1,AKT2,AKT3,ALK,APC,AR,ARAF,AREG,ARID1A,ARID2,ASXL1,
             ATM,ATR,ATRX,AURKA,B2M,BAP1,BCL2,BRAF,BRCA1,BRCA2,BTK,CALR,
             CBL,CCND1,CCND2,CCND3,CCNE1,CD274,CD79B,CDH1,CDK4,CDK6,CDKN2A,CDKN2B,
             CHEK2,CREBBP,CRKL,CSF1R,CTNNB1,CXCR4,DDIT3,DDR2,DICER1,DNAJB1,DNMT3A,EGFR,
             EPCAM,ERBB2,ERBB3,ERBB4,EREG,ERG,ERRFI1,ESR1,EZH1,EZH2,FANCA,FBXW7,
             FGF19,FGF3,FGF4,FGFR1,FGFR2,FGFR3,FGFR4,FLT3,FOXL2,GATA2,GATA3,GNA11,
             GNAQ,GNAS,H3F3A,H3F3B,HGF,HNF1A,HRAS,HSD3B1,IDH1,IDH2,IGF1R,IGF2,
             IKZF1,JAK2,JAK3,KCNJ5,KDM6A,KDR,KEAP1,KIT,KMT2C,KRAS,MAP2K1,MAP2K2,
             MAPK1,MDM2,MED12,MET,MLH1,MSH2,MSH6,MTOR,MUTYH,MYC,MYCN,MYD88,
             NF1,NF2,NFE2L2,NOTCH1,NPM1,NRAS,NTRK1,NTRK2,NTRK3,PALB2,PAX5,PBRM1,
             PDCD1LG2,PDGFRA,PDGFRB,PIK3CA,PIK3CB,PIK3CD,PIK3R1,PML,PMS2,POLD1,POLE,PPP2R1A,
             PPP2R2A,PRKACA,PRKD1,PTCH1,PTEN,PTPN11,RAC1,RAD51B,RAD51C,RAD51D,RAF1,RB1,
             RET,RHEB,RHOA,RIT1,RNF43,ROS1,RUNX1,SERPINB3,SERPINB4,SF3B1,SMAD4,SMARCA4,
             SMARCB1,SMO,SPOP,SRSF2,STAG2,STAT3,STK11,TERT,TET2,TGFBR2,TOP2A,TP53,
             TSC1,TSC2,U2AF1,VEGFA,VEGFB,VHL,WT1,ZNF148,
             CYP2C19,CYP2D6,CYP3A4,DPYD,GSTP1,MTHFR,UGT1A1,XPC,XRCC1,
             AKT2,AREG,AURKA,BCL2,CCNE1,CDKN2B,EREG,FGF19,FGF3,FGF4,HGF,IGF1R,
             IGF2,TOP2A,PPP2R2A,
             ALK,CD74,DDIT3,EGFR,FGFR1,FGFR2,FGFR3,FUS,MET,NTRK1,PDGFRA,PDGFRB,
             PML,PRKACA,RET,ROS1,STIL)))[-1]

genelist_rna <- Hmisc::Cs(
  ABL1,AKT1,ALK,AMACR,ARHGAP26,AXL,BAIAP2L1,BCL2,BCL6,BRAF,CCND1,
  CREB3L2,CREBBP,CTNNB1,DDIT3,DEK,DUSP22,EGFR,ELK4,EPOR,ERBB2,ERG,ESR1,
  ESRRA,EWSR1,FGFR1,FGFR2,FGFR3,FOXO1,GLI1,GLIS2,KLF2,KMT2A,LPP,MALT1,
  MAML2,MAST2,MET,MKL1,MLF1,MYB,MYC,MYH11,NCOA2,NF1,NFIB,NFKB2,
  NOTCH2,NR4A3,NRG1,NTRK1,NTRK2,NTRK3,NUP214,P2RY8,PAX3,PAX5,PAX8,PBX1,
  PDCD1LG2,PDGFB,PDGFRA,PDGFRB,PKN1,PLAG1,PPARG,PRKACA,PRKCA,RAF1,RARA,RELA,
  RET,ROS1,RSPO2,RUNX1,RUNX1T1,SS18,SSX1,SSX2,STAT6,STIL,SUZ12,TAL1,
  TFE3,TFEB,TP63,WDFY2,
  ABL1,ABL2,ABL3,ABL4,ABL5,ABL6,ABL7,ABL8,ABL9,ABL10,ARID1A,ARID2,
  ASXL1,ATM,ATR,AURKA,AXIN2,BAP1,BCL2,BCLAF1,BIRC5,BMPR2,BRAF,BRCA1,
  BRCA2,BTK,CBL,CCND1,CCND2,CCND3,CCNE1,CD274,CD79B,CD8A,CDK6,CDKN2A,
  CDKN2B,CEBPA,CHEK2,CREBBP,CRKL,CTCF,CTNNB1,DICER1,DNMT3A,DOCK3,EGFR,EP300,
  EPCAM,ERBB2,ERBB3,EREG,ESR1,ESRP1,EZH1,EZH2,FBXW7,FGF19,FGF3,FGF4,
  FGFR1,FGFR2,FGFR3,FGFR4,FLT3,FOXL2,GATA2,GATA3,GLI1,GNA11,GNAQ,GPRIN2,
  H3F3A,HNF1A,HRAS,IDH1,IDH2,IGF1R,JAK2,JAK3,KCNJ12,KCNJ5,KIT,KLF4,
  KMT2C,KRAS,MADCAM1,MAP2K1,MAPK1,MED12,MET,MLH1,MLL2,MLL3,MPL,MSH3,
  MSH6,MTOR,MUTYH,MYC,MYCN,MYD88,NEFH,NF2,NFE2L2,NFKBIE,NOTCH1,NOTCH2,
  NR3C1,NRAS,NTRK1,PALB2,PARP4,PAX5,PDCD1,PDCD1LG2,PDGFRA,PDGFRB,PGM5,PIK3CA,
  PIK3R1,PLCG1,PMS2,POLD1,POLE,PPP2R2A,PRIM2,PRKACA,PRKCB,PTCH1,PTEN,PTPN11,
  RAB7A,RAC1,RAD51B,RAD51C,RAD51D,RAF1,RB1,RET,RHEB,RIT1,RNF43,ROS1,
  RUNX1,SERPINB4,SETBP1,SLFN11,SMAD4,SMARCA4,SMO,SPOP,SRSF2,STAG2,STAT3,TET2,
  TGFBR2,TMPRSS13,TOP2A,TP53,TSC1,TSC2,TSHR,U2AF1,UBR5,USP6,VEGFA,VEGFB,
  VHL,WT1,XPO1,XYLT2
)

pathview(gene.data = genelist,
         gene.idtype = 'SYMBOL'
         )



gene.df <- bitr(genelist, fromType = 'SYMBOL', toType = c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db
                )
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                       toType="kegg",organism='hsa')

kegg <- enrichKEGG(gene.df$ENTREZID,
                   organism = 'hsa',
                   keyType = "kegg",
                   pvalueCutoff=0.5,
                   pAdjustMethod='BH',
                   qvalueCutoff=1
                   )
keggx <- setReadable(kegg,'org.Hs.eg.db','ENTREZID')

barplot(kegg, showCategory = 20)
dotplot(kegg, showCategory = 20)

write_delim(keggx@result, 'D:/A549_kegg.tsv', delim = '\t')


# GO
ego <- enrichGO(gene = genelist,
                OrgDb = org.Hs.eg.db,
                ont = 'BP',
                pAdjustMethod='BH', 
                pvalueCutoff=0.05,
                qvalueCutoff=0.2, 
                keyType='SYMBOL'
                )
barplot(ego, showCategory = 20) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )
  
dotplot(ego, showCategory = 20)
egosimp <- simplify(ego,cutoff=0.7,
                    by="p.adjust",
                    select_fun=min,
                    measure="Wang")

write_delim(ego@result, 'D:/A549_GO.tsv', delim = '\t')
write_delim(egosimp@result, 'D:/A549_GO_simply.tsv', delim = '\t')






panel72 <- Hmisc::Cs(ARID1A,
                    HSD3B1,
                    MDM4,
                    AKT3,
                    MSH2,
                    MSH6,
                    ERCC3,
                    NFE2L2,
                    IDH1,
                    FANCD2,
                    MLH1,
                    CTNNB1,
                    FOXP1,
                    RYBP,
                    PIK3CB,
                    ATR,
                    PIK3CA,
                    FBXW7,
                    PIK3R1,
                    CHD1,
                    APC,
                    FANCE,
                    CDK6,
                    MET,
                    BRAF,
                    CUL1,
                    KMT2C,
                    NKX3-1,
                    CLU,
                    NCOA2,
                    MYC,
                    CDKN2A,
                    FANCG,
                    FANCC,
                    PTEN,
                    FANCF,
                    CCND1,
                    ATM,
                    ZBTB16,
                    CDKN1B,
                    KRAS,
                    KMT2D,
                    CDK4,
                    MDM2,
                    BRCA2,
                    RB1,
                    ERCC5,
                    FOXA1,
                    RAD51B,
                    AKT1,
                    IDH2,
                    ERCC4,
                    ZFHX3,
                    FANCA,
                    TP53,
                    CDK12,
                    BRCA1,
                    SPOP,
                    RNF43,
                    RAD51C,
                    AKT2,
                    ERCC2,
                    ERCC1,
                    ASXL1,
                    GNAS,
                    RUNX1,
                    ERG,
                    TMPRSS2,
                    KDM6A,
                    AR,
                    MED12,
                    SMARCA1)




# qianxiang ---------------------------------------------------------------

# dep <- readxl::read_excel('D:/qianxiang0607/TFRE-Kintor-Myc-Protac.xlsx',
#                           sheet = 'c vs d_Myc_1.5nM'
#                           )
# foo <- dep %>% filter(!is.na(DEPs))
tfre <- readxl::read_excel('D:/qianxiang0607/MYC_TFRE.xlsx', 
                           sheet = 'ttest') %>% 
  mutate(deg3nm = case_when(
    ratio_3nM_ctrl > 1.5 & ttest_3nM_ctrl <= 0.05 ~ 'up',
    ratio_3nM_ctrl < 0.67 & ttest_3nM_ctrl <= 0.05 ~ 'down',
    TRUE ~ 'N'
  ), deg15nm = case_when(
    ratio_1.5nM_ctrl > 1.5 & ratio_1.5nM_ctrl <= 0.05 ~ 'up',
    ratio_1.5nM_ctrl < 0.67 & ratio_1.5nM_ctrl <= 0.05 ~ 'down',
    TRUE ~ 'N'
  ))
foo <- tfre %>% 
  filter(ratio_3nM_ctrl > 1.5 | ratio_3nM_ctrl < 0.67) %>% 
  filter(ttest_3nM_ctrl <= 0.05)
# 1.5nM
foo <- tfre %>% 
  filter(ratio_1.5nM_ctrl > 1.5 | ratio_1.5nM_ctrl < 0.67) %>% 
  filter(ttest_1.5nM_ctrl <= 0.05)

read.gmt <- function(path){
  myc <- read_delim(path,
                    delim = '\t',
                    col_names = F
  ) %>% dplyr::select(-X2) %>% 
    column_to_rownames('X1') %>% 
    t() %>% as.data.frame()
  myc
}

myc_v1 <- read.gmt('D:/qianxiang0607/HALLMARK_MYC_TARGETS_V1.gmt')
myc_v2 <- read.gmt('D:/qianxiang0607/HALLMARK_MYC_TARGETS_V2.gmt')
myc_q2 <- read.gmt('D:/qianxiang0607/MYC_Q2.gmt')
myc_h.all <-read.gmt('D:/qianxiang0607/h.all.v7.3.symbols.gmt')

myc_targets <- read_delim('D:/qianxiang0607/MYC_targets.human.tsv', 
                          delim = '\t')
myc_regulators <- read_delim('D:/qianxiang0607/MYC_regulators.human.tsv', 
                          delim = '\t')

myc_ttrust <- base::union(myc_targets$Target, myc_regulators$`# [TF]`)
myc_hallmark <- base::union(myc_v1$HALLMARK_MYC_TARGETS_V1, 
                            myc_v2$HALLMARK_MYC_TARGETS_V2)
myc_transfac <- read_delim('D:/qianxiang0607/MYC_target.txt', 
                           delim = '\t')
  
# 以TTRUST数据库为???
cohort <- '1.5nMTRANSFAC'
degs <- foo$GeneSymbol
degs <- intersect(myc_transfac$gene, degs)
# degs <- intersect(degs, myc_hallmark)
print(length(degs))
write_delim(as.data.frame(degs),
            file.path('D:/qianxiang0607',paste0(cohort,'_DEGs.txt')),
            delim = '\t'
            )

gene.df <- bitr(degs, fromType = 'SYMBOL', 
                toType = c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db
)

kegg <- enrichKEGG(gene.df$ENTREZID,
                   organism='hsa',
                   keyType="ncbi-geneid",
                   pvalueCutoff=0.05,
                   pAdjustMethod='BH',
                   qvalueCutoff=0.2,
                   use_internal_data=F
)
keggx <-  setReadable(kegg,'org.Hs.eg.db',
                      'ENTREZID')

kegg@pvalueCutoff <- 1
kegg@qvalueCutoff <- 1

kegg_p <- barplot(kegg, showCategory = 20,
                  title = 'KEGG Analysis'
                  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )
dotplot(kegg, showCategory = 20)

ggsave(file.path('D:/qianxiang0607',paste0(cohort,'_kegg.pdf')),
       kegg_p
       )
write_delim(keggx@result, 
            file.path('D:/qianxiang0607',paste0(cohort,'_kegg.txt')),
            delim = '\t')

# GO
ego <- enrichGO(gene = degs,
                OrgDb = org.Hs.eg.db,
                ont = 'BP',
                pAdjustMethod='BH', 
                pvalueCutoff=0.05,
                qvalueCutoff=0.2, 
                keyType='SYMBOL'
)

ego@pvalueCutoff <- 1
ego@qvalueCutoff <- 1

go_p <- dotplot(ego, showCategory = 20,
                title = 'Enrichment GO(BP)'
                ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))

egosimp <- simplify(ego,
                    cutoff=0.7,
                    by="p.adjust",
                    select_fun=min,
                    measure="Wang")

ggsave(file.path('D:/qianxiang0607',paste0(cohort,'_go.pdf')),
       go_p)
write_delim(ego@result, 
            file.path('D:/qianxiang0607',paste0(cohort,'_GO.txt')), 
            delim = '\t')
write_delim(egosimp@result, 
            file.path('D:/qianxiang0607',paste0(cohort,'_GO_Sim.txt')),
            delim = '\t')

#DO富集分析
library(DOSE)
do <- enrichDO(gene = gene.df$ENTREZID,
               ont = "DO",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 500,
               qvalueCutoff = 0.2,
               readable = TRUE)

barplot(do, showCategory=20,
        title="Enrichment DO")

#Ractome
x <- enrichPathway(gene=gene.df$ENTREZID,
                   organism = "human",
                   pvalueCutoff=0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.2, 
                   readable=FALSE
                   )

################### string ########################
library(STRINGdb)
library(igraph)
library(ggraph)

string_db <- STRINGdb$new(species=9606,
                          version='11',
                          score_threshold=400
                          )
STRINGdb$methods()
data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "SYMBOL", 
                                      removeUnmappedRows = TRUE)
print(dim(data_mapped))
string_db$plot_network(data_mapped$STRING_id)
data_links <- data_mapped$STRING_id %>% string_db$get_interactions()

links <- data_links %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)

nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% 
  distinct()
net <- igraph::graph_from_data_frame(d=links,
                                     vertices=nodes,
                                     directed = F)

igraph::V(net)$deg <- igraph::degree(net)
igraph::V(net)$size <- igraph::degree(net)/5
igraph::E(net)$width <- igraph::E(net)$weight/10
# 使用ggraph绘图
ggraph(net,layout = "stress")+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>0, label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()

# 如果links数据框的一个link的from只出现过一次，同时to也只出现一次，则将其去除
links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, 
                                                            count(., from)$from)]) %>%
  mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
  filter(!(from_c == 1 & to_c == 1)) %>%
  dplyr::select(1,2,3)
# 新的节点数据
nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络???
net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
# 添加必要的参???
igraph::V(net_2)$deg <- igraph::degree(net_2)
igraph::V(net_2)$size <- igraph::degree(net_2)/5
igraph::E(net_2)$width <- igraph::E(net_2)$weight/10

ppi_p <- ggraph(net_2,layout = "centrality", cent = deg)+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>0, label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()

ggsave(file.path('D:/qianxiang0607',paste0(cohort,'_PPI.pdf')),
       ppi_p,
       width = 15,
       height = 15
       )






