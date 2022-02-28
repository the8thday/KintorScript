

library(tidyverse)
library(ComplexHeatmap)
library(corrplot)
library(magrittr)



# cell_anno <- read_delim('D:/ccle_work/CCLE_sample_info_file_2012-10-18.txt', 
#                         delim = '\t'
#                         )
cell_anno <- read_delim('D:/ccle_work/Cell_lines_annotations_20181226.txt',
                        delim = '\t'
                        )
dpmapID <- cell_anno %>% select(CCLE_ID, depMapID)

ccle_fusion <- read_delim('D:/ccle_work/CCLE_Fusions_20181130.txt',
                          delim = '\t'
                          )
ccle_depmap_maf <- read_delim('D:/ccle_work/CCLE_DepMap_18q3_maf_20180718.txt',
                              delim = '\t'
                              )


depmap_ccle_expression <- vroom::vroom('D:/ccle_work/CCLE_expression.csv') #depmap用的TPM
# depmap_ccle_mutations <- read_csv('D:/ccle_work/CCLE_mutations.csv')
# depmap_ccle_expression %>% filter(X1=='ACH-000440') %>% select(starts_with('MYC'))
depmap_cell_anno <- read_csv('D:/ccle_work/Depmap_sample_info.csv')

# rsem_gene_tpm %>% separate(col = gene_id, into = c('gene_id',NA)) %>% 
#   filter(gene_id=='ENSG00000136997') %>% 
#   select(starts_with('CA46'))


genelist <- Hmisc::Cs(TP53, MYC)
gene.df <- clusterProfiler::bitr(genelist, fromType = 'SYMBOL', 
                                 toType = c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db::org.Hs.eg.db
)


# begin analysis ----------------------------------------------------------

batch1 <- read_delim('D:/ccle_work/batch1.txt',
                     delim = '\t'
                     )
batch2 <- read_delim('D:/ccle_work/batch2.txt',
                     delim = '\t'
)
batch3 <- read_delim('D:/ccle_work/batch3.txt',
                     delim = '\t'
)
batch4 <- read_delim('D:/ccle_work/batch4.txt',
                     delim = '\t'
                     )

drug <- c('STSP', 'GT19715', 'GT19716', 'GT19718')

# 3批次所用细胞系不同
batch_all <- bind_rows(batch1, batch2, batch3, batch4)



f_batch <- function(batch, drug){
  batch %>% filter(`Compound ID`==drug) %>% 
    arrange(`Abs IC50 (uM)`)
}

foo <- f_batch(batch = batch_all, 'GT19715')
psych::describe(foo, omit = T)
Hmisc::describe(foo)

cell_anno_sub <- cell_anno %>% 
  select(CCLE_ID, Name, Histology, Hist_Subtype1, Site_Primary, Disease, type, 
         PATHOLOGIST_ANNOTATION
         )

foo$`Cell Name`[which(foo$`Cell Name`=='Z-138')] <- 'Z138'
foo$`Cell Name`[which(foo$`Cell Name`=='NCI-H526 [H526]')] <- 'NCI-H526'
foo$`Cell Name`[which(foo$`Cell Name`=='NCI-H146 [H146]')] <- 'NCI-H146'
foo$`Cell Name`[which(foo$`Cell Name`=='NCI-H889 [H889]')] <- 'NCI-H889'
foo$`Cell Name`[which(foo$`Cell Name`=='MM.1S')] <- 'MM1-S'
foo$`Cell Name`[which(foo$`Cell Name`=='Jiyoye')] <- 'JiyoyeP-2003'
foo$`Cell Name`[which(foo$`Cell Name`=='SCLC.21H')] <- 'SCLC-21H'
foo$`Cell Name`[which(foo$`Cell Name`=='H2141')] <- 'NCI-H2141'
foo$`Cell Name`[which(foo$`Cell Name`=='H69')] <- 'NCI-H69'
foo$`Cell Name`[which(foo$`Cell Name`=='H9')] <- 'NCI-H9'

jj <- foo %>% left_join(cell_anno_sub,
                 by = c('Cell Name'='Name')
                 )


# 按照细胞系分成两大类
jj %<>% mutate(class = case_when(
  Histology == 'carcinoma' ~ 'SCC',
  TRUE ~ 'Lym'
)) %>% 
  separate(col = CCLE_ID, into = 'ccle', sep = '_') %>% 
  mutate(ccle = if_else(is.na(ccle), `Cell Name`, ccle))
jj$`Cell Name` <- factor(jj$`Cell Name`, levels = jj$`Cell Name`)

write_delim(jj, 'D:/ccle_work/GT19718_anno.txt', delim = '\t')

jj %>% filter(class == 'SCC') %>% 
  ggplot(., aes(x = `Cell Name`, y = `Abs IC50 (uM)`)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))


# analysis by DepMap ------------------------------------------------------

depmap_cell_anno_sub <- depmap_cell_anno %>% select(
  DepMap_ID, cell_line_name, stripped_cell_line_name, primary_disease, Subtype,
  lineage, lineage_subtype,sample_collection_site
)

jj <- foo %>% left_join(depmap_cell_anno_sub,
                        by = c('Cell Name'='cell_line_name')
)

write_delim(jj, 'D:/ccle_work/GT19718_anno_DepMap.txt', delim = '\t')



# mutation ----------------------------------------------------------------
# 所用文件为CCLE_DepMap_18q3_maf_20180718.txt

# input cell lines interesting
intrest_cell <- Hmisc::Cs(EB2, DAUDI) %>% str_to_upper() # top
intrest_cell <- Hmisc::Cs(HEL, Raji, IM9) %>% str_to_upper() # tail
# 只针对small_cell_carcinoma
intrest_cell <- Hmisc::Cs(NCIH2081, DMS53) %>% str_to_upper() #并无交集

ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, '^RAJI')) %>% 
  select(Hugo_Symbol,Protein_Change, cDNA_Change, Tumor_Sample_Barcode) %>% 
  unite(col = 'key', c(Hugo_Symbol, cDNA_Change))


inter_cells <- function(cells){
  stopifnot(length(cells)>=2)
  g1 <- ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, cells[1])) %>% 
    select(Hugo_Symbol,Protein_Change, cDNA_Change) %>% 
    unite(col = 'key', c(Hugo_Symbol, cDNA_Change))
  g1 <- g1$key
  for (i in 2:length(cells)) {
    g2 <- ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, cells[i])) %>% 
      select(Hugo_Symbol,Protein_Change, cDNA_Change) %>% 
      unite(col = 'key', c(Hugo_Symbol, cDNA_Change))
    g1 <- intersect(g1, g2$key)
  }
  return(g1)
}

inter_cells(intrest_cell)


# use depmap to analysis
depmap_ccle_mutations <- vroom::vroom('D:/ccle_work/CCLE_mutations.csv')

depmap_id <- c('ACH-000763', 'ACH-001064') # MM1S, EB2
inter_cells <- function(cells){
  stopifnot(length(cells)>=2)
  g1 <- depmap_ccle_mutations %>% filter(DepMap_ID==cells[1]) %>% 
    select(Hugo_Symbol,Protein_Change ) %>% 
    unite(col = 'key', Hugo_Symbol:Protein_Change)
  g1 <- g1$key
  for (i in 2:length(cells)) {
    g2 <- depmap_ccle_mutations %>% filter(DepMap_ID==cells[i]) %>% 
      select(Hugo_Symbol,Protein_Change ) %>% 
      unite(col = 'key', Hugo_Symbol:Protein_Change)
    g1 <- intersect(g1, g2$key)
  }
  return(g1)
}

inter_cells(intrest_cell)

# 找到所有34个细胞系中的MYC突变
aa <- ccle_depmap_maf %>% separate(col = Tumor_Sample_Barcode, into = 'ccle') %>% 
  filter(Hugo_Symbol=='MYC') %>% 
  filter(ccle %in% c(jj$ccle, c('CTB1','JIYOYEP2003','MAVER1','NCIH9','IM9'))) %>% 
  select(Hugo_Symbol, ccle, Variant_Classification, cDNA_Change, Protein_Change)
write_delim(aa, 'D:/ccle_work/drug1.txt', delim = '\t')

############ 找不敏感细胞系比之其他所特有的位点

i_cells <- Hmisc::Cs(Z138,
                     MM1S,
                     EB2,
                     SUPT1,
                     DAUDI,
                     U937,
                     GA10,
                     C8166,
                     NAMALWA,
                     CTB1,
                     GRANTA519,
                     MM.1R,
                     Ramos,
                     JIYOYEP2003,
                     RL,
                     JEKO1,
                     CA46,
                     PFEIFFER,
                     MAVER1,
                     KARPAS422) %>% str_to_upper()
i_cells <- Hmisc::Cs(NCIH2029,NCIH1092,NCIH526,NCIH82,NCIH146) %>% str_to_upper()

union_cells <- function(cells){
  stopifnot(length(cells)>=2)
  g1 <- ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, cells[1])) %>% 
    select(Hugo_Symbol,Protein_Change, cDNA_Change) %>% 
    unite(col = 'key', c(Hugo_Symbol, cDNA_Change))
  g1 <- g1$key
  for (i in 2:length(cells)) {
    g2 <- ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, cells[i])) %>% 
      select(Hugo_Symbol,Protein_Change, cDNA_Change) %>% 
      unite(col = 'key', c(Hugo_Symbol, cDNA_Change))
    g1 <- union(g1, g2$key)
  }
  return(g1)
}

aa <- ccle_depmap_maf %>% filter(str_detect(Tumor_Sample_Barcode, '^NCIH2081')) %>% 
  select(Hugo_Symbol,Protein_Change, cDNA_Change, Tumor_Sample_Barcode) %>% 
  unite(col = 'key', c(Hugo_Symbol, cDNA_Change))

top_sensitive <- union_cells(i_cells)

setdiff(aa$key,unique(top_sensitive))

library(VennDiagram)
venn_list <- list(other=unique(top_sensitive), NCIH2081=aa$key)
venn.plot <- venn.diagram(
  venn_list,
  filename = 'venn.png',
  fill = c('red', 'green'),
  alpha = 0.5,
  col = 'transparent',
  print.mode = c('row'),
  output = TRUE
)


# FUSION ------------------------------------------------------------------

interest_cells <- Hmisc::Cs(NCIH2081,DMS53)

a2 <- ccle_fusion %>% filter(str_detect(`X.sample`, '^RAJI')) %>% 
  select(`X.sample`,`X.FusionName`,LeftBreakpoint,RightBreakpoint) %>% 
  distinct()

find_fusion_overlap <- function(cells){
  stopifnot(length(cells)>=2)
  g1 <- ccle_fusion %>% filter(str_detect(`X.sample`, cells[1])) %>% 
    select(`X.sample`,`X.FusionName`,LeftBreakpoint,RightBreakpoint) %>% 
    distinct()
  g1 <- g1$`X.FusionName`
  for (i in 2:length(cells)) {
    g2 <- ccle_fusion %>% filter(str_detect(`X.sample`, cells[i])) %>% 
      select(`X.sample`,`X.FusionName`,LeftBreakpoint,RightBreakpoint) %>% 
      distinct()
    g1 <- intersect(g1, g2$X.FusionName)
  }
  return(g1)
}

find_fusion_overlap(interest_cells)



# methylation -------------------------------------------------------------

ccle_methy <- vroom::vroom('D:/ccle_work/CCLE_RRBS_TSS1kb_20181022.txt/CCLE_RRBS_TSS1kb_20181022.txt')
# ccle_methy_2 <- read_delim('D:/ccle_work/CCLE_RRBS_TSS1kb_20181022.txt/CCLE_RRBS_TSS1kb_20181022.txt',
#                            delim = '\t'
#                            )

ccle_methy %>% select(starts_with('MM1S'), locus_id)

find_methy_overlap <- function(cells){
  g1 <- ccle_methy %>% select(starts_with(cells[1]), locus_id)
  g1 <- g1$locus_id
  for (i in 2:length(cells)) {
    g2 <- ccle_methy %>% select(starts_with(cells[i]), locus_id)
    g1 <- intersect(g1, g2$locus_id)
  }
  print(length(g1))
}

interest_cells <- Hmisc::Cs(NCIH2081, DMS53)

find_methy_overlap(interest_cells)


# 对于甲基化应该求其相关性

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

ccle_methy %>% select(starts_with('NCIH2081') | starts_with('DMS53'), locus_id) %>%
  drop_na() %>% 
  column_to_rownames('locus_id') %>% cor() %>% round(., 3) %>% 
  get_lower_tri(.) %>% 
  reshape2::melt(data = ., na.rm = TRUE) %>% 
  ggplot(., aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = 'white')+
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "white",
                       midpoint = 0.5, limit = c(0.1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


# RNA expression ----------------------------------------------------------

rpkm <- read_delim('D:/ccle_work/CCLE_RNAseq_genes_rpkm_20180929.gct/CCLE_RNAseq_genes_rpkm_20180929.gct',
                   delim = '\t',
                   skip = 2
)

# 挑选某基因在细胞系中的表达
rpkm %>% filter(Description=='GAPDH') %>% 
  select(starts_with(c('CALU','HAE')))

# 目前来看rpkm算是和网站结果比较一致的
# rsem_gene_tpm <-  read_delim('D:/ccle_work/CCLE_RNAseq_rsem_genes_tpm_20180929.txt/CCLE_RNAseq_rsem_genes_tpm_20180929.txt',
#                              delim = '\t'
#                              )
rna_broad <- vroom::vroom('D:/ccle_work/rnase')

genelist <- c('TP53','MYC','BCL2','BCL6','KRAS','JAK1','JAK2','FGFR1','FGFR2',
              'FGFR3','RUNX1','CREBBP','EP300','MALT1','BCR'
              )
# genelist_ren <- unique(Hmisc::Cs(NRAS,KRAS,FLT3,IL7R,SH2B3,BRAF,GATA3,ETV6,RUNX1, EP300,PAX5,RB1,JAK2,CDKN2B,CDKN2A,NOTCH2,
#                           TP53,SF3B1,ATM,BIRC3,NOTCH1,KLF2,BCL6,TNFAIP3,CD79A,EZH2,ARID1A,MEF2B,PTEN,GNA13,B2M,CD58,
#                           CREBBP,EP300,FOXO1,KMT2C,CCND1,CARD11,PTPRD, HER2, IGFR, EGFR, FGFR,PI3K,STAT3,STAT5B,ATM,
#                           JAK1,MAPK1,CXCR4,TCF3,ID3,BCL6,MAP2K1,NOTCH2,KMT2D,KLF2,SPEN, IDH1,IDH2, CD79B))
# genelist <- genelist_ren

allcells <- jj %>% filter(class == 'SCC')
allcells <- allcells$ccle %>% str_to_upper() %>% 
  str_c('_')

M <- rpkm %>% 
  # filter(Description %in% genelist) %>% 
  select(Description, starts_with(allcells)) %>%
  # select(matches(paste(allcells, collapse="|")))
  distinct(Description, .keep_all = T) %>%
  column_to_rownames('Description')

# Z SCORE
cn <- colnames(M) %>% map(., function(x)str_split(x, '_')[[1]][1]) %>% 
  as.vector()
colnames(M) <- cn

Heatmap(log2(M + 1),
        column_names_gp = grid::gpar(fontsize = 12),
        column_names_rot = 60,
        cluster_columns = F,
        cluster_rows = T,
        row_names_gp = grid::gpar(fontsize = 6)
        )
cor(log2(M + 1))


# depmap
depmap_ccle_expression %>% filter(`...1` %in% c(''))

###################### 求每个基因的表达值和药物IC50 的相关性 #############

drug <- jj %>% filter(class == 'SCC') %>%
  select(ccle, `Abs IC50 (uM)`) %>% 
  filter(ccle %in% colnames(M)) %>% 
  rename('IC50'='Abs IC50 (uM)') %>% 
  column_to_rownames('ccle') %>% 
  t() %>% 
  as.data.frame()

cor.test(as.numeric(M['BCL6',]), as.numeric(drug['IC50',]))

outTab <- data.frame()
M %<>% drop_na()
M <- M[rowSums(M==0)<5,]

for(Gene in row.names(M)){
  x <- as.numeric(M[Gene,])
  for(Drug in row.names(drug)){
    y <- as.numeric(drug[Drug,])
    corT <- cor.test(x,y,method="pearson")
    cor <- corT$estimate
    pvalue <- corT$p.value
    if(pvalue < 0.01){
      outVector <- cbind(Gene,Drug,cor,pvalue)
      outTab <- rbind(outTab,outVector)
    }
  }
}

outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]


library(ggpubr)

plotList_1 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}

for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(M[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4])<0.001){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1=ggplot(data = df1, aes(x = x, y = y))+
    geom_point(size=1)+
    stat_smooth(method="lm",se=FALSE, formula=y~x)+
    labs(x="Expression",y="IC50",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
    theme_bw()
  plotList_1[[i]]=p1
}

plotList_2 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}


for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2]
  x <- as.numeric(M[Gene,])
  y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(x,y))
  colnames(df1)[2] <- "IC50"
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  compaired <- list(c("low", "high"))
  p1 <- ggboxplot(df1,
                  x = "group", y = "IC50",
                  fill = "group", palette = c("#00AFBB", "#E7B800"),
                  add = "jitter", size = 0.5,
                  xlab = paste0("Expression_of_", Gene),
                  ylab = paste0("IC50_", Drug)) +
    stat_compare_means(comparisons = compaired,
                       method = "wilcox.test",   #设置统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))
  plotList_2[[i]]=p1
}



nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)
box <- ggarrange(plotlist=plotList_2,nrow=nrow,ncol=ncol)

ggsave(plot = box, filename = 'D:/ccle_work/box.pdf',
       height = 15,
       width = 15
       )



# drug sensitive ----------------------------------------------------------

drug <- read_delim('D:/ccle_work/primary-screen-replicate-collapsed-logfold-change.csv',
                   delim = ','
                   )
drug1 <- readxl::read_excel('D:/ccle_work/GDSC1_fitted_dose_response_25Feb20.xlsx')
drug2 <- readxl::read_excel('D:/ccle_work/GDSC2_fitted_dose_response_25Feb20.xlsx')

cell34 <- c('Z138','MM1S','EB2','NCI-H526','Daudi','U-937','NCI-H82','NCI-H146',	'C8166',	
            'NAMALWA',	'CTB-1','GRANTA-519','NCI-H841','NCI-H446','MM-1R','JiyoyeP-2003','DMS-273','RL',
                    'JEKO-1','CA46',	'SCLC-21H','Pfeiffer','MAVER-1','NCI-H889',	'NCI-H2171',
                    'KARPAS-422','KMS-11','NCI-H1694','MC116','HEL','IM-9','DMS-53','NCI-H2081','Raji')

i_cell <- c('Raji', 'HEL', 'REC-1', 'IM-9')
i_cell <- c('MM1S','EB2','CTB-1')
i_cell <- c('DMS-53','NCI-H2081')
i_cell <- c('NCI-H526','NCI-H82', 'NCI-H146', 'NCI-H841','NCIH2029', 'NCI-H1092')
i_cell <- c('Raji', 'GRANTA-519','NCI-H9','KMS-11')
i_cell <- c('NCI-H524','NCI-H526')
i_cell <- c('NCI-H446','COR-L47','NCI-H2141','COR-L88','HCI-H69','SCLC-21H')

bar <- drug2 %>% filter(CELL_LINE_NAME %in% i_cell) %>% 
  arrange(Z_SCORE) %>% 
  select(CELL_LINE_NAME, TCGA_DESC, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, Z_SCORE) %>% 
  filter(Z_SCORE < -2) %>% arrange(PUTATIVE_TARGET)

write_delim(bar, 'D:/ccle_work/drug1.txt', delim = '\t')
  


# CNV ---------------------------------------------------------------------

# CNV 所使用文件为
cnv <- read_csv('D:/ccle_work/cnv_20191101/cnv_abs_copy_number_picnic_20191101.csv',
                skip = 1
                )

cnv %>% select(`MM1S`, `Raji`)
# PICNIC 为SNP6.0芯片数据分析



ccle_translocations <- readxl::read_excel('D:/ccle_work/CCLE_translocations_SvABA_20181221.xlsx')
# 对于转坐是只关注MYC，还是寻求某种差异？
ccle_translocations %>%
  filter(str_detect(CCLE_name, paste(allcells, collapse="|")))
# 
# depmap_cnv <- vroom::vroom('D:/ccle_work/CCLE_gene_cn.csv')

ccle_translocations %>% filter(gene1=='MYC'|gene2=='MYC')

