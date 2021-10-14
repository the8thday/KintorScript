library(tidyverse)
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


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
    dplyr::select(geneID, pval, FDR, logFC, Regulation, GeneSymbol)
}


R1R2_R2R3 <- go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R3_sig), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  select(sort(names(.)))

R1R2_R2R3 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/R1R2_R2R3.txt',
                          delim = '\t')


R1R2_R2R4 <- go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R4_sig), by = 'geneID', suffix = c('_R1R2', '_R2R4')) %>%
  select(sort(names(.)))
  # write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/R1R2_R2R4.txt', delim = '\t')

R2R3_R2R4 <- go_inter(R2R3_sig) %>%
  full_join(go_inter(R2R4_sig), by = 'geneID', suffix = c('_R2R3', '_R2R4')) %>%
  dplyr::select(sort(names(.)))
  # write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/R2R3_R2R4.txt', delim = '\t')

# merge all three files
foo <- go_inter(R2R4_sig) %>%
  rename_with(.fn = ~paste0(.x, '_R2R4')) %>%
  rename(geneID=geneID_R2R4)
go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R3_sig), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  full_join(foo, by = 'geneID') %>%
  select(sort(names(.)))

# plot heatmap
library(ComplexHeatmap)

M <- R1R2_R2R3 %>%
  select(GeneSymbol_R1R2, GeneSymbol_R2R3, logFC_R1R2, logFC_R2R3) %>%
  filter(!is.na(GeneSymbol_R2R3)) %>%
  arrange(desc(logFC_R1R2)) %>%
  select(-GeneSymbol_R1R2) %>%
  column_to_rownames('GeneSymbol_R2R3')

Heatmap(M,
        name = 'logFC',
        # na_col = '',
        border_gp = gpar(col = 'black'),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        # top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        # column_names_gp = gpar(fontsize=1)
)


## -----U937--------------

path <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Human/Report/Result/08_DE'

U5U6_sig <- read_delim(file.path(path, 'Group_U-5-VS-U-6_DE_significant_anno.xls'),
                       delim = '\t')
U6U7_sig <- read_delim(file.path(path, 'Group_U-6-VS-U-7_DE_significant_anno.xls'),
                       delim = '\t')
U6U8_sig <- read_delim(file.path(path, 'Group_U-6-VS-U-8_DE_significant_anno.xls'),
                       delim = '\t')

go_inter <- function(x){
  x %>%
    dplyr::select(geneID, pvalue, qvalue, log2FoldChange, Regulation, GeneSymbol)
}

U5U6_U6U7 <- go_inter(U5U6_sig) %>%
  full_join(go_inter(U6U7_sig), by = 'geneID', suffix = c('_U5U6', '_U6U7')) %>%
  dplyr::select(sort(names(.)))

U5U6_U6U7 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/U5U6_U6U7.txt',
                          delim = '\t')

U5U6_U6U8 <- go_inter(U5U6_sig) %>%
  full_join(go_inter(U6U8_sig), by = 'geneID', suffix = c('_U5U6', '_U6U8')) %>%
  dplyr::select(sort(names(.)))

U5U6_U6U8 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/U5U6_U6U8.txt',
                          delim = '\t')




# DIY ----------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

genelist <- R1R2_sig %>%
  dplyr::select(geneID, GeneSymbol) %>%
  pull(geneID)

gene.df <- bitr(genelist,
                fromType = 'ENSEMBL',
                toType = c('ENTREZID', 'SYMBOL'),
                OrgDb = org.Mm.eg.db,
                drop = TRUE
)

# go 富集分析
enrich.go <- enrichGO(
  gene = genelist,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 1,
  qvalueCutoff  = 1
)
# enrich.go@pvalueCutoff <- 1
# enrich.go@qvalueCutoff <- 1

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
              title = 'GO(BP) Analysis'
) +
  # scale_color_gradient2() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

pp2 <- barplot(enrich.go,
               font.size = 12
) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

(pp3 <- cnetplot(enrich.go,
                 node_label = 'all'))

ego2 <- enrichplot::pairwise_termsim(enrich.go)
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
  gene = gene.df$ENTREZID,
  keyType = 'ncbi-geneid',
  organism = 'hsa',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  use_internal_data = T
)
kegg@pvalueCutoff <- 1
kegg@qvalueCutoff <- 1

keggx <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')

kegg_p <- barplot(kegg,
                  font.size = 15,
                  title = 'KEGG Analysis'
)
write_delim(kegg@result,
            file = file.path(outpath, paste0(cohort, '_KEGG.txt')), delim = '\t')
write_delim(keggx@result,
            file = file.path(outpath, paste0(cohort, '_KEGGX.txt')), delim = '\t')
ggsave(file.path(outpath, paste0(cohort, '_KEGG.pdf')),
       plot = kegg_p, width = 10, height = 10)

# browseKEGG(kegg, pathID = 'hsa04080')

# WikiPathways, Reactome
ewiki <- enrichWP(gene.df$ENTREZID,
                  organism = "Homo sapiens",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5
)
ewiki_p <- dotplot(ewiki,
                   font.size = 12,
                   color = 'p.adjust',
                   showCategory=20,
                   title = 'WikiPathway ORA Analysis'
) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_WikiPathway.pdf')),
       plot = ewiki_p, width = 10, height = 10)
ewiki <- setReadable(ewiki, OrgDb = org.Hs.eg.db)
write_delim(ewiki@result,
            file = file.path(outpath, paste0(cohort, '_WikiPathway.txt')), delim = '\t')

eReactome <- ReactomePA::enrichPathway(gene=gene.df$ENTREZID,
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
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

# viewPathway("E2F mediated regulation of DNA replication",
#             readable = TRUE,
#             foldChange = geneList)
ggsave(file.path(outpath, paste0(cohort, '_Reactome.pdf')),
       plot = react_p, width = 10, height = 10)
write_delim(eReactome@result,
            file = file.path(outpath, paste0(cohort, '_Reactome.txt')), delim = '\t')


y <- gsePathway(gene.df$ENTREZID,
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH",
                verbose = FALSE)

# DOSE
enrich.do <- DOSE::enrichDO(
  gene = gene.df$ENTREZID,
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
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_DO.pdf')),
       plot = do_p, width = 10, height = 10)
write_delim(enrich.do@result,
            file = file.path(outpath, paste0(cohort, '_DO.txt')), delim = '\t')


# gsea 富集 免疫通路
get_id <- function(dataset){
  dataset %>%
    left_join(uniprt_id, by = c('uniprot_id'='Entry')) %>%
    dplyr::select(score, Symbol) %>%
    left_join(gene.df, by = c('Symbol'='SYMBOL')) %>%
    dplyr::select(ENTREZID, score) %>%
    arrange(desc(score)) %>% distinct()
}
ss <- get_id(gcnnet %>% dplyr::filter(score>-0.97)) # 注意修改数据集
gene_sort <- ss %>% pull(score)
names(gene_sort) <- ss$ENTREZID

c7 <- msigdbr::msigdbr(species = "Homo sapiens", category = 'C7') %>%
  dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
colnames(c7)<-c("ont","gene")

gsea <- GSEA(gene_sort,
             TERM2GENE = c7,
             pvalueCutoff = 0.5)
em2 <- setReadable(gsea,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")
em2.df <- as.data.frame(em2)
em2.df <- em2.df %>%
  dplyr::select(ID,NES,setSize,p.adjust)
GSEA_result2 <- merge(c7,em2.df, by.x='ont', by.y = 'ID')
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

# gsea kegg
gseKEGG.res <- gseKEGG(gene_sort,
                       organism = "hsa",
                       keyType = "kegg",
                       minGSSize = 5,
                       maxGSSize = 500,
                       pvalueCutoff=1
)
gseKEGG_p <- enrichplot::gseaplot2(gseKEGG.res ,1:5,
                                   pvalue_table = FALSE
)

ggsave(file = file.path(outpath, paste0(cohort, '_gseKEGG.pdf')),
       plot = gseKEGG_p, width = 10, height = 10)
write_delim(gseKEGG.res@result,
            file = file.path(outpath, paste0(cohort, '_gseKEGG.txt')), delim = '\t')






# GO KEGG --------------------------------------------------------------------

path2 <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Mouse/Report/Result/09_GO'
path3 <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Mouse/Report/Result/10_KEGG'

R1R2_all_go <- readxl::read_excel(file.path(path3, 'R-1-VS-R-2.xlsx'))
R1R2_all_go <- readxl::read_excel(file.path(path3, 'R-2-VS-R-3.xlsx'))
R1R2_all_go <- readxl::read_excel(file.path(path3, 'R-2-VS-R-4.xlsx'))

koid <- c('ko05171', 'ko05321', 'ko04933', 'ko04940',
          'ko04060','ko04630','ko04668','ko04061'
          )
koname <- c('Coronavirus disease-COVID-19', 'Inflammatory bowel disease',
            'AGE-RAGE signaling pathway in diabetic complications',
            'Type I diabetes mellitus',
            'Cytokine-cytokine receptor interaction',
            'JAK-STAT signaling pathway',
            'TNF signaling pathway',
            'Viral protein interaction with cytokine and cytokine receptor'
            )
names(koid) <- koname

find_gene <- function(dir, ko){
  path4 <- file.path(path3, dir)
  ff <- list.files(
    path = path4,
    pattern = paste0('*.',ko, '_[A-Za-z]{3,4}.xlsx'),
    recursive = TRUE,
    ignore.case = TRUE
  )
  if(length(ff)==0){print(glue::glue('{dir} Has No Such PathWay!'))}
  down <- readxl::read_excel(file.path(path4, ff[1])) %>%
    mutate(Chr=as.character(Chr))
  up <- readxl::read_excel(file.path(path4, ff[2])) %>%
    mutate(Chr=as.character(Chr))
  if(dim(down)==0 && dim(up)!=0){
    return(up %>% mutate(Compare = {{dir}}))
  } else if(dim(down)!=0 && dim(up)==0){
    return(down %>% mutate(Compare = {{dir}}))
  } else {
    return(bind_rows(down, up) %>%
             mutate(Compare = {{dir}}))
  }
}


bar <- koid[2]
find_gene('R-1-VS-R-2', bar) %>%
  bind_rows(find_gene('R-2-VS-R-3', bar)) %>%
  bind_rows(find_gene('R-2-VS-R-4', bar)) %>%
  relocate(Compare, everything()) %>%
  write_delim(file = file.path('/Users/congliu/OneDrive/kintor/Daily_Work',
                               paste0(names(bar),'.txt')),
              delim = '\t')


for(i in seq_along(koid)){
  print(koid[i])
  find_gene('R-1-VS-R-2', koid[i]) %>%
    bind_rows(find_gene('R-2-VS-R-3', koid[i])) %>%
    bind_rows(find_gene('R-2-VS-R-4', koid[i])) %>%
    relocate(Compare, everything()) %>%
    write_delim(file = file.path('/Users/congliu/OneDrive/kintor/Daily_Work',
                                 paste0(names(koid[i]),'.txt')),
                delim = '\t')
}


# another file format
find_gene2 <- function(dir, ko){
  find_gene(dir = dir, ko = ko) %>%
    dplyr::select(Compare, geneID, logFC, pval, FDR, Regulation, GeneSymbol)
}

foo <- find_gene2('R-2-VS-R-4', koid[1]) %>%
  rename_with(.fn = ~paste0(.x, '_R2R4')) %>%
  dplyr::rename('geneID'='geneID_R2R4')

find_gene2('R-1-VS-R-2', koid[1]) %>%
  full_join(find_gene2('R-2-VS-R-3', koid[1]), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  full_join(foo, by = 'geneID') %>%
  dplyr::select(sort(names(.))) %>%
  write_delim(file = file.path('/Users/congliu/OneDrive/kintor/Daily_Work',
                               paste0(names(bar),'_merge.txt')),
              delim = '\t')








