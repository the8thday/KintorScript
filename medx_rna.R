## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2021-11-22
##
## Copyright (c) cliu, 2021
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 这是一批次毛囊RNASEQ的测序数据，由迈杰（medx）提供.
## 共25例样本，其中group12为10对配对样本，group13为5例正常vs10例枕部患者
##
## ---------------------------
library(tidyverse)


path <- '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/result/4_differential_expression'
setwd(path)

all_de <- list.files(
  path = path,
  pattern = '*Genes_DE.txt',
  recursive = TRUE
)



# filter DEGs -------------------------------------------------------------

filter_de <- function(f){
  sample_name <- str_extract(f, 'group[0-9]{1,2}')
  print(sample_name)
  ff <- read_delim(f, delim = '\t',show_col_types = FALSE)
  ff %>% filter(significant=='yes') %>%
    select(gene_name,gene_id,log2FoldChange,pvalue,padj,regulation,significant) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    write_csv(
      glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/de_filter/{sample_name}_filter.csv'
                 )
    )
}

filter_de_edgeR <- function(f){
  sample_name <- str_extract(f, 'group[0-9]{1,2}')
  print(sample_name)
  ff <- read_delim(f, delim = '\t', show_col_types = FALSE)
  ff %>% filter(significant=='yes') %>%
    select(gene_name,gene_id,logFC,PValue,FDR,regulation,significant) %>%
    arrange(desc(abs(logFC))) %>%
    write_csv(
      glue::glue('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/de_filter/{sample_name}_filter.csv')
    )
}

for(i in all_de){
  if(str_detect(i, 'group1/|group12|group13')){
    filter_de(i)
  } else {
    filter_de_edgeR(i)
  }
}


# gsea --------------------------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Pi)
library(viridis)

outpath = '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/de_filter/'
cohort = 'Group12'

df_de <- read_delim(
  file = 'group12/group12Genes_DE.txt',
  delim = '\t'
)

ss <- df_de %>%
  # filter(significant=='yes') %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  group_by(gene_name) %>%
  top_n(1, abs(log2FoldChange)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(desc(log2FoldChange))

gene_map <- clusterProfiler::bitr(ss$gene_name,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
ss <- ss %>% left_join(gene_map, by = c('gene_name'='SYMBOL')) %>%
  drop_na()

gene_sort <- ss %>% pull(log2FoldChange)
names(gene_sort) <- ss$ENTREZID
gene_sort <- na.omit(gene_sort)


# gmt gsea
methy <- msigdbr::msigdbr(species = "Homo sapiens",
                       category = 'C5')
methy_term <- methy %>%
  # filter(gs_name == 'GOBP_DNA_METHYLATION') %>%
  filter(str_detect(gs_name,'METHYLATION')) %>%
  dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

methy_id <- methy %>%
  filter(str_detect(gs_name,'METHYLATION')) %>%
  dplyr::select(gs_name, entrez_gene, gene_symbol, ensembl_gene, gs_description)
df_de %>% filter(gene_name %in% unique(methy_id$gene_symbol)) %>%
  write_delim(file = file.path(outpath, paste0(cohort, '_methy.txt')),
              delim = '\t')

gsea <- GSEA(gene_sort,
             TERM2GENE = methy_term,
             pvalueCutoff = 0.2
             )
em2 <- setReadable(gsea,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")
em2.df <- as.data.frame(em2)
em2.df <- em2.df %>%
  dplyr::select(ID,NES,setSize,p.adjust)
GSEA_result <- merge(c7,em2.df, by.x='ont', by.y = 'ID')
id <- 1
anno <- em2[id, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

gsea_p1 <- enrichplot::gseaplot(gsea, 1,
                                title = gsea$Description[1]
) +
  annotate("text", 0.8, 0.85,
           label = lab, hjust=0, vjust=0)

gsea_p <- enrichplot::gseaplot2(gsea, 1:2,
                                title = 'GSEA Analysis'
)
ridgeplot(gsea)
ggsave(file = file.path(outpath, paste0(cohort, '_GSEA.pdf')),
       plot = gsea_p, width = 10, height = 10)


# gsea go
gseGo.res <- gseGO(geneList=gene_sort,
             ont ="BP",
             keyType = "ENTREZID",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             by = 'fgsea',
             nPermSimple = 10000
             )

gseGo.res <- simplify(gseGo.res,
                      cutoff = 0.7,
                      by = "p.adjust",
                      select_fun = min,
                      measure = "Wang")

gseGO_p <- enrichplot::gseaplot2(gseGo.res,
                      1:5,
          title = "GSEA GO",
          color="red",
          base_size = 20,
          subplots = 1:3,
          pvalue_table = T
          )
ggsave(file = file.path(outpath, paste0(cohort, '_gseGO.pdf')),
       plot = gseGO_p, width = 20, height = 10)

gseGo.resX <- setReadable(gseGo.res,
                         OrgDb = org.Hs.eg.db)
write_delim(gseGo.resX@result,
            file = file.path(outpath, paste0(cohort, '_gseGO.txt')), delim = '\t')

gseGO_p2 <- dotplot(gseGo.res,
        showCategory = 20,
        title = 'GSEA GO Analysis'
        ) +
  scale_color_viridis()
  # scale_color_continuous(
  #   low="red", high="blue",
  #   guide=guide_colorbar(reverse=TRUE)
  # )

ggsave(file = file.path(outpath, paste0(cohort, '_gseGO_dotplot.pdf')),
       plot = gseGO_p2, width = 10, height = 10)

ggplot(gseGo.res@result %>% head(),
       aes(x=NES,y=Description,color=p.adjust,
                        size=setSize))+
  geom_point()+
  scale_size_continuous(breaks = c(40,80,120,180),
                        name = "No. of\nsignificant genes")+
  scale_color_gradientn(colours = rainbow(4),limits=c(0,1),
                        breaks=c(0,0.25,0.5,0.75,1),
                        labels=c(0,0.25,0.5,0.75,1))+  xlim(c(-3,3))+
  geom_vline(xintercept=0)+  xlab("Normalized enrichment score")+
  ylab("")+  theme_bw()

# gsea kegg
gseKEGG.res <- gseKEGG(gene_sort,
                       organism = "hsa",
                       keyType = "kegg",
                       minGSSize = 5,
                       maxGSSize = 500,
                       pvalueCutoff=0.05,
                       pAdjustMethod = "BH"
)
gseKEGG_p <- enrichplot::gseaplot2(gseKEGG.res, 1:5,
                                   pvalue_table = FALSE
)

ggsave(file = file.path(outpath, paste0(cohort, '_gseKEGG.pdf')),
       plot = gseKEGG_p, width = 10, height = 10)

gseKEGG.resX <- setReadable(gseKEGG.res,
                            keyType = 'ENTREZID',
                            OrgDb = org.Hs.eg.db)
write_delim(gseKEGG.resX@result,
            file = file.path(outpath, paste0(cohort, '_gseKEGG.txt')), delim = '\t')

(gseKEGG_p2 <- dotplot(gseKEGG.res,
                    showCategory = 20,
                    title = 'GSEA KEGG Analysis'
) +
  scale_color_viridis())

ggsave(file = file.path(outpath, paste0(cohort, '_gseKEGG_dotplot.pdf')),
       plot = gseKEGG_p2, width = 10, height = 10)


# Pi xPierGSEA, xPierPathways
ss_pi <- ss %>% pull(log2FoldChange)
names(ss_pi) <- ss$gene_name
ss_pi <- na.omit(ss_pi)

Pnode <- xPierGenes(ss_pi)


pnode_ss <- ss %>% filter(!gene_name %in% c('MEMO1','TEC')) %>%
  column_to_rownames('gene_name')
names(pnode_ss) <- c('rank', 'priority')

eGSEA <- Pi::xPierGSEA(
  pNode = pnode_ss,
  ontology = "MsigdbC2KEGG",
  size.range = c(10, 500),
  nperm = 10000,
  fast = FALSE,
  p.adjust.method = 'BH',
  RData.location = "http://galahad.well.ox.ac.uk/bigdata"
  # RData.location = '/Users/congliu/OneDrive/kintor/xPierGSEA_download/ontology_Rdata'
  )

xGSEAdotplot(eGSEA, top=1)


# KEGG GO -----------------------------------------------------------------

# 差异基因GO，KEGG分析
de_list <- df_de %>%
  filter(significant=='yes') %>%
  dplyr::select(gene_name, gene_id)
gene_map <- clusterProfiler::bitr(de_list$gene_id,
                                  fromType = "ENSEMBL",
                                  toType = c("ENTREZID", 'SYMBOL', 'REFSEQ'),
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
de_l <- unique(gene_map$ENTREZID)


enrich.go <- enrichGO(
  gene = de_list$gene_name,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.1,
  qvalueCutoff  = 0.2,
  # readable = TRUE
)
enrich.go <- setReadable(enrich.go,
                          OrgDb = org.Hs.eg.db)

egosim <- simplify(enrich.go,
                   cutoff = 0.7,
                   by = 'p.adjust',
                   select_fun = min,
                   measure = 'Wang'
)

pp <- dotplot(enrich.go,
              font.size = 12,
              color = 'p.adjust',
              showCategory=20,
              title = 'GO(BP) Enrichment'
) +
  scale_color_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

pp2 <- barplot(egosim,
               font.size = 12
) +
  scale_fill_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

(pp3 <- cnetplot(egosim,
                 node_label = 'all'))

ego2 <- enrichplot::pairwise_termsim(egosim)
(pp4 <- enrichplot::treeplot(ego2,
                             hclust_method = "ward.D"
))

write_delim(enrich.go@result,
            file = file.path(outpath, paste0(cohort, '_goResult.txt')),
            delim = '\t')
write_delim(egosim@result,
            file = file.path(outpath, paste0(cohort, '_goSim.txt')),
            delim = '\t')
ggsave(file = file.path(outpath, paste0(cohort, '_goResult.pdf')),
       plot = pp, width = 10, height = 10)
ggsave(file = file.path(outpath, paste0(cohort, '_goCnet.pdf')),
       plot = pp3, width = 10, height = 10)
ggsave(file = file.path(outpath, paste0(cohort, '_goTree.pdf')),
       plot = pp4, width = 10, height = 10)

# kegg
kegg <- enrichKEGG(
  gene = de_l,
  keyType = 'ncbi-geneid',
  organism = 'hsa',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  use_internal_data = F
)
kegg@pvalueCutoff <- 1
kegg@qvalueCutoff <- 1

keggx <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')

kegg_p <- barplot(kegg,
                  font.size = 15,
                  title = 'KEGG Enrichment',
                  showCategory = 20
) +
  scale_fill_viridis()

write_delim(keggx@result,
            file = file.path(outpath, paste0(cohort, '_KEGGX.txt')), delim = '\t')
ggsave(file.path(outpath, paste0(cohort, '_KEGG.pdf')),
       plot = kegg_p, width = 10, height = 10)

# browseKEGG(kegg, pathID = 'hsa04080')

# WikiPathways, Reactome
ewiki <- enrichWP(de_l,
                  organism = "Homo sapiens",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2
)
ewiki_p <- dotplot(ewiki,
                   font.size = 12,
                   color = 'p.adjust',
                   showCategory=20,
                   title = 'WikiPathway ORA Enrichment'
) +
  scale_color_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_WikiPathway.pdf')),
       plot = ewiki_p, width = 10, height = 10)
ewiki <- setReadable(ewiki, OrgDb = org.Hs.eg.db)
write_delim(ewiki@result,
            file = file.path(outpath, paste0(cohort, '_WikiPathway.txt')), delim = '\t')

eReactome <- ReactomePA::enrichPathway(gene=de_l,
                                       pvalueCutoff = 0.05,
                                       readable=TRUE,
                                       organism = "human",
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.2
)
react_p <- dotplot(eReactome,
                   font.size = 12,
                   color = 'p.adjust',
                   showCategory=20,
                   title = 'Reactome Analysis'
) +
  scale_color_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

# viewPathway("E2F mediated regulation of DNA replication",
#             readable = TRUE,
#             foldChange = geneList)
ggsave(file.path(outpath, paste0(cohort, '_Reactome.pdf')),
       plot = react_p, width = 10, height = 10)
write_delim(eReactome@result,
            file = file.path(outpath, paste0(cohort, '_Reactome.txt')), delim = '\t')


y <- gsePathway(de_l,
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH",
                verbose = FALSE)

# DOSE
enrich.do <- DOSE::enrichDO(
  gene = de_l,
  ont = 'DO',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

do_p <- dotplot(enrich.do,
                font.size = 12,
                color = 'p.adjust',
                showCategory=20,
                title = 'DO Analysis'
) +
  scale_color_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_DO.pdf')),
       plot = do_p, width = 10, height = 10)
write_delim(enrich.do@result,
            file = file.path(outpath, paste0(cohort, '_DO.txt')), delim = '\t')


# WGCNA -------------------------------------------------------------------

library(WGCNA)







# AR gene & follicle geneset GSEA analysis---------------------------------------

library(magrittr)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))


path1 <- '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA'

all_fpkm <- read_delim(
  '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/result/3_expression/Gene_fpkm_expression.txt',
  delim = '\t'
  )
all_fpkm %<>% drop_na(gene_id)
all_fpkm <- all_fpkm %>% rowwise() %>%
  mutate(s = sum(c_across(starts_with('FPKM')))) %>%
  filter(s > 0) %>% dplyr::select(-s)
sample_info <- readxl::read_excel(file.path(path1,'sample25.xlsx'))

AR_Gene1 <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/AR_targets.human.tsv') %>%
  dplyr::rename('TF'='# TF')
AR_Gene2 <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/AR_regulators.human.tsv')
names(AR_Gene2) <- names(AR_Gene1)
AR_Gene <- AR_Gene1 %>% bind_rows(AR_Gene2)

ar_go <- read.gmt(file.path(path1, 'GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
wp_go <- read.gmt(file.path(path1, 'WP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
ar_reg <- read.gmt(file.path(path1, 'GOBP_REGULATION_OF_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
ar_msigdb <- ar_go %>% bind_rows(wp_go) %>% bind_rows(ar_reg)
# ar_id <- unique(ar_msigdb %>% pull(gene))

ar_id <- unique(c(unique(AR_Gene %>% pull(TF)), unique(AR_Gene$Target)))
ar_id <- unique(c(AR_Gene2$TF, 'AR'))
ar_id <- unique(hair_one$gene)
print(length(ar_id))
rm_normal <- sample_info %>% filter(class !='normal') %>%
  mutate(bianhao = str_c('FPKM.',bianhao)) %>%
  pull(bianhao)

M <- all_fpkm %>% dplyr::filter(`gene_name` %in% ar_id) %>%
  dplyr::select(`gene_name`, starts_with('FPKM')) %>%
  # dplyr::select(`gene_name`, rm_normal) %>%
  # filter(gene_name != 'PIK3R2') %>%
  # mutate(`FPKM.12506N190026`=as.numeric(`FPKM.12506N190026`)) %>%
  column_to_rownames('gene_name')
skimr::skim(M)

s_info <- sample_info %>% dplyr::select(-sample) %>%
  # filter(class !='normal') %>%
  column_to_rownames('bianhao')
ha <- HeatmapAnnotation(df = s_info, col =
                          list(df = c("normal" = "red", "occ" = "green", "v" = "blue")))

ComplexHeatmap::Heatmap(log2(M+1),
        name = 'ratio',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = T,
        show_column_dend = T,
        cluster_columns = T,
        cluster_rows = T,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        column_names_gp = gpar(fontsize=8),
        row_names_gp = gpar(fontsize=8)
)


# GSEA aalysis ---
cohort <- 'group12'
path2 <- '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/'

df_de <- read_delim(
  file = '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/group_12_paired/DESeq2_DEG_lfcShrink_group12.txt',
  delim = '\t'
)

ss <- df_de %>%
  # filter(significant=='yes') %>%
  dplyr::select(ENTREZID, log2FoldChange) %>%
  drop_na() %>%
  group_by(ENTREZID) %>%
  top_n(1, abs(log2FoldChange)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(desc(log2FoldChange))

gene_map <- clusterProfiler::bitr(ss$ENTREZID,
                                  fromType = "ENTREZID",
                                  toType = "SYMBOL",
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
ss <- ss %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  left_join(gene_map, by = c('ENTREZID'='ENTREZID')) %>%
  drop_na()

gene_sort <- ss %>% pull(log2FoldChange)
names(gene_sort) <- ss$ENTREZID
gene_sort <- na.omit(gene_sort)
print(length(gene_sort))

AnnotationDbi::select(x = org.Hs.eg.db, keys = 'GO:0051799',
                      columns = c('ENSEMBL', 'GENENAME'),
                      keytype = 'GOALL'
                      )


# AR geneset
ar_go <- read.gmt(file.path(path1, 'GOBP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
wp_go <- read.gmt(file.path(path1, 'WP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
ar_reg <- read.gmt(file.path(path1, 'GOBP_REGULATION_OF_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY.gmt'))
ar_ttrust <- read.gmt(file.path(path1, 'AR_TTRUST_GENE.gmt'))
ar_msigdb <- ar_go %>% bind_rows(wp_go) %>% bind_rows(ar_reg) %>% bind_rows(ar_ttrust)
ar_one <- ar_msigdb %>%
  mutate(term = 'MSigDB_TRRUST_AR') %>% distinct()

# 筛选AR基因相关的差异基因
target_file <- df_de %>% dplyr::filter(SYMBOL %in% unique(ar_msigdb$gene))
write_delim(x = target_file,
            file = file.path(path2, paste0(cohort, '_AR.txt')),
            delim = '\t'
)

# hair gmt gsea
go_mature <- read.gmt(file.path(path1, 'GOBP_HAIR_FOLLICLE_MATURATION.gmt'))
go_neg <- read.gmt(file.path(path1, 'GOBP_NEGATIVE_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT.gmt'))
go_dev <- read.gmt(file.path(path1, 'GOBP_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT.gmt'))
go_pos <- read.gmt(file.path(path1, 'GOBP_POSITIVE_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT.gmt'))
wp_dev <- read.gmt(file.path(path1, 'WP_HAIR_FOLLICLE_DEVELOPMENT_CYTODIFFERENTIATION_PART_3_OF_3.gmt'))
wp_dev2 <- read.gmt(file.path(path1, 'WP_HAIR_FOLLICLE_DEVELOPMENT_ORGANOGENESIS_PART_2_OF_3.gmt'))
wp_eda <- read.gmt(file.path(path1, 'WP_EDA_SIGNALLING_IN_HAIR_FOLLICLE_DEVELOPMENT.gmt'))
GOBP_METHYLATION <- read.gmt('/Users/congliu/Downloads/geneset.gmt')

hair_all <- go_mature %>%
  bind_rows(go_neg) %>%
  bind_rows(go_dev) %>%
  bind_rows(wp_dev2) %>%
  bind_rows(go_pos) %>%
  bind_rows(wp_dev) %>%
  bind_rows(wp_eda)

hair_one <- hair_all %>%
  filter(term != 'GOBP_NEGATIVE_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT') %>%
  mutate(term = 'MSigDB_HAIR_FOLLICLE') %>% distinct()

change_gmt <- function(x){
  gene_map2 <- clusterProfiler::bitr(x$gene,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
  x <- x %>% left_join(gene_map2, by = c('gene'='SYMBOL')) %>%
    dplyr::select(term, ENTREZID) %>% dplyr::rename('gene'='ENTREZID') %>%
    drop_na()
  x
}

geneset <- change_gmt(hair_all)
universe <- gene_map$ENTREZID

gsea <- GSEA(gene_sort,
             TERM2GENE = geneset,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             minGSSize = 3,
             maxGSSize = 500,
             by = "fgsea"
)
em2 <- setReadable(gsea,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")

em2.df <- as.data.frame(em2)
GSEA_result <- merge(geneset,em2.df, by.x='term', by.y = 'ID')

id <- 1
anno <- em2[id, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

gsea_p1 <- enrichplot::gseaplot(em2, 1,
                                title = gsea$Description[1]
) +
  annotate("text", 0.8, 0.85,
           label = lab, hjust=0, vjust=0)

gsea_p <- enrichplot::gseaplot2(gsea, 1,
                                title = gsea$Description[1]
) +
  annotate("text", 0.8, 0.85,
           label = lab, hjust=0, vjust=0)

# ridgeplot(gsea)
ggsave(file = file.path(path2, paste0(cohort, '_GSEA.pdf')),
       plot = gsea_p1, width = 10, height = 10)
ggsave(file = file.path(path2, paste0(cohort, '_GSEA2.pdf')),
       plot = gsea_p, width = 10, height = 10)

write_delim(x = em2@result,
            file = file.path(path2, paste0(cohort, '_GSEA.txt')),
            delim = '\t'
            )

leading_gene <- (em2@result %>% pull(core_enrichment) %>%
  str_split(pattern = '/'))[[1]]
gsea_de <- df_de %>% filter(SYMBOL %in% leading_gene)
  dplyr::select(SYMBOL,log2FoldChange, pvalue, padj, regulation, significant)

write_delim(x = gsea_de,
            file = file.path(path2, paste0(cohort, '_GSEA_DE.txt')),
            delim = '\t'
)

id <- em2@geneSets
geneInCategory(em2)[id]


## enricher
ss2 <- df_de %>%
  filter(significant=='yes') %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  group_by(gene_name) %>%
  top_n(1, abs(log2FoldChange)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(desc(log2FoldChange))

gene_map <- clusterProfiler::bitr(ss2$gene_name,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
ss2 <- ss2 %>% left_join(gene_map, by = c('gene_name'='SYMBOL')) %>%
  drop_na()

gene_de <- ss2 %>% pull(ENTREZID)

em <- enricher(gene_de,
               TERM2GENE=geneset,
               pvalueCutoff = 0.5,
               pAdjustMethod = 'BH',
               qvalueCutoff = 1,
               universe = universe
               )
em <- setReadable(em,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")
head(em)


# plot DEGs by logFC ------------------------------------------------------

library(ComplexHeatmap)
library("RColorBrewer")

target_genes <- read_delim(
  '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/group_12_paired/AR/group12_AR.txt'
) %>%
  drop_na()
# target_genes <- target_genes %>% filter(SYMBOL %in% c(unique(AR_Gene1$Target), 'AR'))
target_genes <- target_genes %>% filter(SYMBOL %in% c(unique((target_genes %>% slice(1:33))$SYMBOL), 'AR'))

fpkm <- read_tsv('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/FPKM.txt')
metafile <- "/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/sample25.xlsx"

fpkm <- fpkm %>%
  # dplyr::rename('rowname'='Gene ID') %>%
  mutate(rowname = gsub("\\.\\d*", "",rowname))
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm <- apply(fpkm %>% column_to_rownames('rowname'),2,fpkmToTpm)
# fpkm <- as.data.frame(tpm) %>%
#   rownames_to_column()

fpkm <- fpkm %>%
  rowwise() %>% dplyr::mutate(sum = sum(c_across(where(is.numeric)))) %>%
  filter(sum>0) %>% dplyr::select(-sum)
sampleGroup <- readxl::read_excel(metafile) %>%
  arrange(class)


# select id
ids_remove <- sampleGroup %>% filter(class=='occ') %>% pull(bianhao)
ids_remove <- c(ids_remove, '12506N0001'
                # '12506I190015','12506N0002','12506I190014'
)

fpkm <- fpkm %>%
  filter(rowname %in% target_genes$ENSEMBL) %>%
  dplyr::select(!{{ids_remove}}) %>%
  distinct()

anno1 <- sampleGroup %>% dplyr::select(bianhao, class, sample) %>%
  filter(class!='occ') %>%
  filter(!bianhao %in% c('12506N0001'))

foo <- fpkm %>% column_to_rownames('rowname') %>% t() %>% as.data.frame() %>%
  rownames_to_column() %>% left_join(anno1, by = c('rowname'='bianhao')) %>%
  arrange(class) %>% dplyr::select(-class, -rowname) %>% rename('rowname'='sample') %>%
  column_to_rownames('rowname') %>% t() %>% as.data.frame() %>%
  rownames_to_column() %>% left_join(target_genes,by = c('rowname'='ENSEMBL')) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(-c(rowname,ENTREZID,baseMean,lfcSE,stat,pvalue,padj))

anno <- columnAnnotation(group = anno1$class,
                         col = list(group = c('normal'='#009999', 'v'='coral1'))
                         )

ha <- rowAnnotation(logFC = anno_barplot(foo$log2FoldChange,
                                         gp = gpar(fill = 'darksalmon',
                                                   col = 'darksalmon'
                                                   )))

M <- foo %>% dplyr::select(-log2FoldChange) %>%
  column_to_rownames('SYMBOL')

col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

png(filename = '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/group13/AR/AR_FPKM_group13.png',
    width=10,height=8,units="in",res=1000
)
ht <- Heatmap(t(scale(t(log2(M + 1)))),
        name = 'FPKM',
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
        row_names_side = "left",
        col = circlize::colorRamp2(c(-4, 0, 4), c("navy", "white", "firebrick3")),
        # col = circlize::colorRamp2(c(-4, 0, 4), brewer.pal(n=3, name="RdBu")),
        width = ncol(M)*unit(5, "mm"),
        # height = nrow(M)*unit(6, "mm"),
        right_annotation = ha,
        top_annotation = anno,
        column_names_rot = 90,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize=8),
        )

dev.off()















