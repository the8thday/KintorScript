## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2021-11-17
##
## Copyright (c) cliu, 2021
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 金唯智RNA测序数据的再分析，为两细胞系U937/RAW，在lps刺激后加入普克鲁胺
##
##
## ---------------------------
library(tidyverse)
library(ComplexHeatmap)
library(VennDiagram)
library(pheatmap)
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


# compare DEGs RAW ------------------------------------------------------------

path <- '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse'
outpath <- path
cohortname <- 'mouse'
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png <- function(x, filename, width=8,height=10,units="in",res=1000) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height,units = units, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

all_fpkm <- read_delim(
  '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/all.fpkm_anno.xls',
  delim = '\t'
  )


R1R2 <- read_delim(file.path(path, 'Group_R-1-VS-R-2_DE_anno.xls'),
                   delim = '\t')
R2R3 <- read_delim(file.path(path, 'Group_R-2-VS-R-3_DE_anno.xls'),
                   delim = '\t')
R2R4 <- read_delim(file.path(path, 'Group_R-2-VS-R-4_DE_anno.xls'),
                   delim = '\t')

R1R2 %>% mutate(GeneSymbol = str_to_upper(GeneSymbol)) %>%
  filter(str_detect(GeneSymbol, paste0('^',i_gene,collapse = '|'))) %>%
  filter(GeneSymbol %in% i_gene) %>%
  select(GeneSymbol, logFC)


R1R2_sig <- read_delim(file.path(path, 'Group_R-1-VS-R-2_DE_significant_anno.xls'),
                       delim = '\t')
R2R3_sig <- read_delim(file.path(path, 'Group_R-2-VS-R-3_DE_significant_anno.xls'),
                       delim = '\t')
R2R4_sig <- read_delim(file.path(path, 'Group_R-2-VS-R-4_DE_significant_anno.xls'),
                       delim = '\t')


go_inter <- function(x){
  x %>%
    dplyr::select(gene_id, PValue, FDR, logFC, GeneSymbol)
}


R1R2_R2R3 <- go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R3_sig), by = 'geneID', suffix = c('_R1R2', '_R2R3')) %>%
  select(sort(names(.)))

R1R2_R2R3 <- go_inter(R1R2) %>%
  full_join(go_inter(R2R3), by = 'gene_id', suffix = c('_R1R2', '_R2R3')) %>%
  select(sort(names(.)))

# R1R2_R2R3 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/R1R2_R2R3.txt',
#                           delim = '\t')


R1R2_R2R4 <- go_inter(R1R2_sig) %>%
  full_join(go_inter(R2R4_sig), by = 'geneID', suffix = c('_R1R2', '_R2R4')) %>%
  select(sort(names(.)))
  # write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/R1R2_R2R4.txt', delim = '\t')

R2R3_R2R4 <- go_inter(R2R3_sig) %>%
  full_join(go_inter(R2R4_sig), by = 'geneID', suffix = c('_R2R3', '_R2R4')) %>%
  dplyr::select(sort(names(.)))
  # write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/R2R3_R2R4.txt', delim = '\t')

# merge all three files
foo <- go_inter(R2R4) %>%
  rename_with(.fn = ~paste0(.x, '_R2R4')) %>%
  rename(gene_id=gene_id_R2R4)
all_three <- go_inter(R1R2) %>%
  full_join(go_inter(R2R3), by = 'gene_id', suffix = c('_R1R2', '_R2R3')) %>%
  full_join(foo, by = 'gene_id') %>%
  select(sort(names(.)))


# plot heatmap

M <- all_three %>%
  select(GeneSymbol_R1R2, GeneSymbol_R2R3, GeneSymbol_R2R4, logFC_R1R2, logFC_R2R3, logFC_R2R4) %>%
  filter(!is.na(GeneSymbol_R2R4)) %>%
  mutate(GeneSymbol_R2R4 = str_to_upper(GeneSymbol_R2R4)) %>%
  filter(GeneSymbol_R2R4 %in% i_gene) %>%
  arrange(desc(logFC_R1R2), desc(logFC_R2R4), desc(logFC_R2R3)) %>%
  # filter(!is.na(GeneSymbol_R1R2)) %>%
  select(-GeneSymbol_R1R2, -GeneSymbol_R2R3) %>%
  filter(!str_detect(GeneSymbol_R2R4, '^Gm')) %>%
  filter(!str_detect(GeneSymbol_R2R4, 'Rik$')) %>%
  dplyr::rename_with(.fn = ~str_extract(.x, 'R[0-9]R[0-9]'), .col=starts_with('logFC')) %>%
  column_to_rownames('GeneSymbol_R2R4')
# M[is.na(M)] <- 0


p_heat1 <- pheatmap::pheatmap(M,
                             cluster_rows = FALSE,
                             show_rownames = TRUE,
                             cluster_cols = FALSE,
                             clustering_method = "average",
                             # scale = "row",
                             cellwidth = 15,
                             # cellheight = 15,
                             fontsize_row = 8,
                             border = FALSE,
                             legend = TRUE,
                             legend_labels = 'logFC'
                             # na_col = 'grey'
                             # annotation_col = anno
)
save_pheatmap_png(p_heat1,
                  file.path(
                    outpath,
                    paste0("heatmap_nodiff_", cohortname, ".png")
                  )
)

save_pheatmap_pdf(p_heat1,
                  file.path(
                    outpath,
                    paste0("heatmap_nodiff_", cohortname, ".pdf")
                  )
                  )

png(file.path(
  outpath,
  paste0("heatmap_nodiff_", cohortname, ".png")),
    width=8,height=10,units="in",res=1000
)
Heatmap(M,
        name = 'logFC',
        col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
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
        column_names_rot = 45,
        show_column_names = T,
        row_names_gp = gpar(fontsize=6),
        width = unit(20, 'mm')
)

dev.off()

# venn plot
g12 <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/paired/DESeq2_DEG_filtered_group12.txt',
                  delim = '\t') %>% drop_na()
g13 <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/group13/DESeq2_DEG_filtered_group13.txt',
                  delim = '\t') %>% drop_na()
venn_list <- list(g12=g12$SYMBOL,g13=g13$SYMBOL)

venn.plot <- venn.diagram(venn_list,
                          filename = '/Users/congliu/OneDrive/kintor/Daily_Work/RAW_venn.png',
                          fill = c('red', 'green'), alpha = 0.50,
                          # cat.col = rep('black', 2),
                          # col = 'black',
                          col = "transparent",
                          print.mode=c("raw","percent"),
                          output=TRUE,
                          # cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif', margin = 0.2,
                          # imagetype="png",
                          # height = 880,
                          # width = 880,
                          # resolution = 1000,
                          # compression = "lzw"
                          )


## -----U937--------------

path <- '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Human'
outpath <- path
cohortname <- 'U937'

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

# U5U6_U6U7 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/U5U6_U6U7.txt',
#                           delim = '\t')

U5U6_U6U8 <- go_inter(U5U6_sig) %>%
  full_join(go_inter(U6U8_sig), by = 'geneID', suffix = c('_U5U6', '_U6U8')) %>%
  dplyr::select(sort(names(.)))
#
# U6U7_U6U8 %>% write_delim(file = '/Users/congliu/OneDrive/kintor/Daily_Work/U6U7_U6U8.txt',
#                           delim = '\t')

U6U7_U6U8 <- go_inter(U6U7_sig) %>%
  full_join(go_inter(U6U8_sig), by = 'geneID', suffix = c('_U6U7', '_U6U8')) %>%
  dplyr::select(sort(names(.)))

# merge all three files
foo <- go_inter(U6U8_sig) %>%
  rename_with(.fn = ~paste0(.x, '_U6U8')) %>%
  dplyr::rename(geneID=geneID_U6U8)
all_three <- go_inter(U5U6_sig) %>%
  full_join(go_inter(U6U7_sig), by = 'geneID', suffix = c('_U5U6', '_U6U7')) %>%
  full_join(foo, by = 'geneID') %>%
  dplyr::select(sort(names(.)))


M <- all_three %>%
  select(GeneSymbol_U5U6, GeneSymbol_U6U7, GeneSymbol_U6U8,log2FoldChange_U5U6, log2FoldChange_U6U7, log2FoldChange_U6U8) %>%
  filter(!is.na(GeneSymbol_U6U8)) %>%
  filter(GeneSymbol_U6U8 != '-') %>%
  mutate(GeneSymbol_U6U8 = str_to_upper(GeneSymbol_U6U8)) %>%
  # filter(GeneSymbol_U6U8 %in% i_gene) %>%
  arrange(desc(log2FoldChange_U5U6), desc(log2FoldChange_U6U8), desc(log2FoldChange_U6U7)) %>%
  filter(GeneSymbol_U6U8 != 'DNAJC9-AS1') %>%
  filter(!is.na(GeneSymbol_U5U6)) %>%
  select(-GeneSymbol_U6U7, -GeneSymbol_U5U6) %>%
  dplyr::rename_with(.fn = ~str_extract(.x, 'U[0-9]U[0-9]'), .col=starts_with('log2FoldChange')) %>%
  column_to_rownames('GeneSymbol_U6U8')

p_heat1 <- pheatmap::pheatmap(M,
                              cluster_rows = FALSE,
                              show_rownames = TRUE,
                              cluster_cols = FALSE,
                              clustering_method = "average",
                              # scale = "row",
                              cellwidth = 15,
                              # cellheight = 15,
                              fontsize_row = 3,
                              border = FALSE,
                              # na_col = 'grey'
                              # annotation_col = anno
)
save_pheatmap_pdf(p_heat1,
                  file.path(
                    outpath,
                    paste0("heatmap_diff_", cohortname, ".pdf")
                  ),
                  width = 10, height = 10
)

pdf("/Users/congliu/OneDrive/kintor/Daily_Work/three_groups_2.pdf",
    width=8,height=8)
Heatmap(M,
        name = 'logFC',
        col = circlize::colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick3")),
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
        column_names_rot = 45,
        show_column_names = T,
        row_names_gp = gpar(fontsize=1),
        width = unit(20, 'mm'),
        column_names_gp = gpar(fontsize=8)
)

dev.off()

# venn plot
venn_list <- list(U5U6_sig=U5U6_sig$geneID,
                  U6U7_sig=U6U7_sig$geneID,
                  U6U8_sig=U6U8_sig$geneID)
venn_list <- list(
  U5U6_sig = U5U6_sig %>% filter(GeneSymbol != '-') %>% pull(geneID),
  U6U7_sig = U6U7_sig %>% filter(GeneSymbol != '-') %>% pull(geneID),
  U6U8_sig = U6U8_sig %>% filter(GeneSymbol != '-') %>% pull(geneID)
)

venn.plot <- venn.diagram(venn_list,
                          filename = '/Users/congliu/OneDrive/kintor/Daily_Work/U937_del_lncRNA_venn.png',
                          fill = c('red', 'green', 'blue'), alpha = 0.50,
                          # cat.col = rep('black', 2),
                          # col = 'black',
                          col = "transparent",
                          print.mode=c("raw","percent"),
                          output=TRUE,
                          # cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif', margin = 0.2,
                          # imagetype="png",
                          # height = 880,
                          # width = 880,
                          # resolution = 1000,
                          # compression = "lzw"
)

# DIY ----------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

cohort <- 'R1R2_R2R4'

# genelist <- R1R2_R2R4 %>%
#   filter(!is.na(GeneSymbol_R2R4)) %>%
#   filter(!is.na(GeneSymbol_R1R2)) %>%
#   dplyr::select(geneID) %>%
#   distinct() %>%
#   pull(geneID)

# genelist <- U5U6_U6U8 %>%
#   filter(!is.na(GeneSymbol_U6U8)) %>%
#   filter(!is.na(GeneSymbol_U5U6)) %>%
#   dplyr::select(geneID) %>%
#   distinct() %>%
#   pull(geneID)

genelist <- R2R4_sig$geneID
geneMap <- clusterProfiler::bitr(genelist,
                                 fromType = "ENSEMBL",
                                 toType = c("SYMBOL", "ENTREZID"),
                                 OrgDb = "org.Mm.eg.db",
                                 drop = TRUE
)


# go 富集分析
enrich.go <- enrichGO(
  gene = genelist,
  OrgDb = org.Mm.eg.db,
  keyType = 'ENSEMBL',
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff  = 0.2
)
# enrich.go@pvalueCutoff <- 1
# enrich.go@qvalueCutoff <- 1

enrich.go <- setReadable(enrich.go,
                         OrgDb = org.Mm.eg.db)
egosim <- simplify(enrich.go,
                   cutoff = 0.7,
                   by = 'p.adjust',
                   select_fun = min,
                   measure = 'Wang'
)

pp <- dotplot(egosim,
              font.size = 12,
              color = 'p.adjust',
              showCategory=20,
              title = glue::glue('{cohortname} GO(BP) Enrichment')
) +
  viridis::scale_color_viridis() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

pp2 <- barplot(egosim,
               font.size = 12
) +
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
  gene = geneMap$ENTREZID,
  keyType = 'ncbi-geneid',
  organism = 'mmu',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  use_internal_data = F
)
# kegg@pvalueCutoff <- 1
# kegg@qvalueCutoff <- 1

keggx <- setReadable(kegg, 'org.Mm.eg.db', 'ENTREZID')

kegg_p <- clusterProfiler::dotplot(kegg,
                                   font.size = 15,
                                   title = glue::glue('{cohortname} KEGG Enrichment'),
                                   showCategory = 20
) +
  scale_color_viridis() +
  scale_y_discrete(labels = function(x) {
    str_wrap(x, width = 40)
  })

kegg_p2 <- barplot(kegg,
                   font.size = 8,
                   showCategory = 20,
                   title = glue::glue('{cohortname} KEGG Enrichment')
) +
  scale_fill_viridis() +
  scale_y_discrete(labels = function(x) {
    str_wrap(x, width = 40)
  })
write_delim(kegg@result,
            file = file.path(outpath, paste0(cohortname, '_KEGG.txt')), delim = '\t')
write_delim(keggx@result,
            file = file.path(outpath, paste0(cohortname, '_KEGGX.txt')), delim = '\t')
ggsave(file.path(outpath, paste0(cohortname, '_KEGG.pdf')),
       plot = kegg_p, width = 10, height = 10)
ggsave(file.path(outpath, paste0(cohortname, '_KEGG_bar.pdf')),
       plot = kegg_p2, width = 10, height = 10)

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






# RAW GO KEGG intersect--------------------------------------------------------------------

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



# U937 GO KEGG intersect------------------------------------------------------------

path2 <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Human/Report/Result/09_GO'
path3 <- '/Volumes/H500G-86/80-762814506/N2111843_80-762814506/Human/Report/Result/10_KEGG'

R1R2_all_go <- readxl::read_excel(file.path(path3, 'U-5-VS-U-6.xlsx'))
R1R2_all_go <- readxl::read_excel(file.path(path3, 'U-6-VS-U-7.xlsx'))
R1R2_all_go <- readxl::read_excel(file.path(path3, 'U-6-VS-U-8.xlsx'))

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

i_ko <- c('ko04080', 'ko04024', 'ko04728')

find_gene('U-6-VS-U-8', 'ko04061') %>%
  write_delim('/Users/congliu/OneDrive/kintor/Daily_Work/ko04061_U6U8.txt', delim = '\t')


find_gene_go <- function(dir, ko){
  path4 <- file.path(path2, dir)
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

i_gene <- c('GO_0006955', 'GO:0051607', 'GO:0006935', 'GO:0009615', 'GO:0009617','GO:0060337',
            'GO:0045071'
            ) %>% str_replace(':', '_')
names(i_gene) <- c('immune response', 'defense response to virus', 'chemotaxis',
                   'response to virus', 'response to bacterium',
                   'type I interferon signaling pathway',
                   'negative regulation of viral genome replication'
                   )

find_gene_go('U-6-VS-U-8', 'GO_0006955')

for(i in seq_along(i_gene)){
  print(i_gene[i])
  find_gene_go('U-6-VS-U-8', i_gene[i]) %>%
    write_delim(file = file.path('/Users/congliu/OneDrive/kintor/Daily_Work',
                                 paste0(names(i_gene[i]),'.txt')),
                delim = '\t')
}















# heatmap by FPKM ---------------------------------------------------------

i_gene <- c(
  "IL1a","IL1RA","IL1b","IL2","IL4","IL5","IL6","IL7","IL9","IL10","IL12b","IL13","IL15","IL17A",
  "IL18","IL21","IL22","IL23a","IL27","IL31","IFNa","IFNg","TNF","CCL2","CCL3","CCL4","CCL5","CCL6",
  "CCL7","CCL22","CCL11","CCL19","CCL20","CCL27","CXCL2","CXCL1","CXCL6","CXCL8","CXCL9","CXCL10",
  "CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17","TNFb","NGFb","BDNF","EGF","FGF2","LIF",
  "PDGFBB","PlGF1","SCF","VEGFA","VEGFD","BAFF","GCSF","GMCSF","MCSF","Csf1r","ICAM1","Csf1","CSF2",
  "CSF3"
) %>% str_to_upper()


all_fpkm <- read_delim(
  '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Human/all.fpkm_anno.xls',
  delim = '\t'
)

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

is.infinite.matrix <- function(x){
  apply(x,2,is.infinite)
}


fpkm_human <- all_fpkm %>%
  select(ends_with('_FPKM'),GeneSymbol) %>%
  rowwise() %>%
  mutate(s = sum(c_across(ends_with('FPKM')))) %>%
  filter(s > 0) %>% dplyr::select(-s) %>%
  mutate(GeneSymbol = str_to_upper(GeneSymbol)) %>%
  filter(GeneSymbol %in% i_gene) %>%
  arrange(desc(`6_FPKM`))
print(dim(fpkm_human))

fpkm_mouse <- all_fpkm %>%
  mutate(GeneSymbol = str_to_upper(GeneSymbol)) %>%
  filter(GeneSymbol %in% i_gene) %>%
  select(ends_with('_FPKM'),GeneSymbol) %>%
  rowwise() %>%
  mutate(s = sum(c_across(ends_with('FPKM')))) %>%
  filter(s > 0) %>% dplyr::select(-s) %>%
  distinct() %>%
  arrange(desc(`1_FPKM`))
print(dim(fpkm_mouse))

tpms <- apply(all_fpkm %>%
                select(ends_with('_FPKM'),gene_id) %>%
                column_to_rownames('gene_id'),
              2,fpkmToTpm)


M <- fpkm_human %>%
  column_to_rownames('GeneSymbol')
print(dim(M))

M <- M[rowSums(M)>1,]

m1 <- log2(M + 1)
m2 <- t(scale(t(M)))

p_heat <- pheatmap::pheatmap(m1,
                   cluster_rows = TRUE,
                   show_rownames = TRUE,
                   cluster_cols = FALSE,
                   clustering_method = "average",
                   scale = "row",
                   cellwidth = 15,
                   cellheight = 15,
                   border = FALSE,
                   # annotation_col = anno
                   )

outpath <- path
cohortname <- 'human'
save_pheatmap_pdf(p_heat,
                  file.path(
                    outpath,
                    paste0("fpkm_", cohortname, ".pdf")
                  ),
                  width = 10, height = 10
)


# ha <- columnAnnotation(bar = anno_barplot(final_data$days))
# pdf()
Heatmap(m1,
        name = 'log10(FPKM+1)',
        na_col = '#E6E6FA',
        border_gp = gpar(col = 'black'),
        show_row_dend = T,
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = T,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        col = circlize::colorRamp2(c(0, 2, 4), c("navy", "white", "firebrick3")),
        # col = circlize::colorRamp2(c(0,16), c("white", "#C00000")),
        width = ncol(M)*unit(8, "mm"),
        # height = nrow(M)*unit(6, "mm"),
        # top_annotation = ha,
        column_names_rot = 90,
        show_column_names = T,
        show_row_names = FALSE,
        # column_names_gp = gpar(fontsize=1)
)

dev.off()













# plot kegg ---------------------------------------------------------------

kegg_df <- read_delim(
  '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/R-2-VS-R-4_DE_significant_enrichment.xls'
)


kegg_df_f <- kegg_df %>% filter(Qvalue<=0.05) %>%
  select(c(Pathway,Level1,`DEGs with pathway annotation (154)`),Qvalue) %>%
  mutate(num = str_extract(`DEGs with pathway annotation (154)`,'[0-9]{1,2}')) %>%
  mutate(num = as.integer(num)) %>%
  arrange(Level1, num)
kegg_df_f$Pathway <- factor(kegg_df_f$Pathway, levels = unique(kegg_df_f$Pathway))

library(showtext)
font_add_google()


(ps <- kegg_df_f %>%
  filter(Level1 %in% c(
                       'Environmental Information Processing'
                       )) %>%
  ggplot(aes(x = Pathway, y = num, fill = Level1)) +
  geom_col(fill = 'firebrick1') +
  geom_text(aes(label=num, y = num+1.5)) +
  # ggrepel::geom_text_repel(aes(label=num)) +
  coord_flip() +
  # ggsci::scale_fill_aaas() +
  theme_bw() +
  scale_x_discrete(labels = function(x) {str_wrap(x, width = 40)}) +
  theme(text = element_text(size = 12),
        legend.position = 'none'
        ) + ggtitle('Environmental Information Processing') )

ggsave(filename = '/Users/congliu/OneDrive/kintor/Daily_Work/genviz_1012/Mouse/ei.pdf',
       plot = ps,
       width = 10, height = 10
         )






