#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

library(tidyverse, quietly = TRUE)
library(DESeq2, quietly = TRUE)
# library(biomaRt, quietly = TRUE)
library("BiocParallel")
# library(Glimma)
register(MulticoreParam(4))
options(stringsAsFactors = F)

metafile <- "/Users/congliu/OneDrive/kintor/qingiqng/GT9001/Group_GT_6h.tsv"
countfile <- "/Users/congliu/OneDrive/kintor/qingiqng/GT9001/all.fpkm_anno.xls"
outpath <- "/Users/congliu/OneDrive/kintor/qingiqng/GT9001/stim_6h/"
cohortname <- "TGF_blank1h"


sampleGroup <- read.table(metafile,
  header = TRUE, row.names = 1
)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("Control", "Test"))

geneCount <- read_delim(countfile, delim = "\t") %>%
  dplyr::select(!ends_with("FPKM")) %>%
  dplyr::select(-`Exonic.gene.sizes`)
geneCount <- geneCount[, -c((dim(geneCount)[2] - 9):(dim(geneCount)[2]))] %>%
  dplyr::select(`gene_id`, any_of(rownames(sampleGroup))) %>%
  column_to_rownames("gene_id")


geneMap <- clusterProfiler::bitr(rownames(geneCount),
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = "org.Hs.eg.db"
)

# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# geneMap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                  filters = "ensembl_gene_id",
#                  values = rownames(geneCount),
#                  mart = ensembl) %>%
#   dplyr::filter(!(is.na(hgnc_symbol) & is.na(entrezgene_id))) %>% dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
print(head(geneMap, n = 3))

# DESeqDataSet/DESeqDataSetFromTximport/DESeqDataSetFromHTSeq
dds <- DESeqDataSetFromMatrix(
  countData = geneCount,
  colData = sampleGroup,
  design = ~Group
)
# 指定哪一组作为control
dds$group <- relevel(dds$Group, ref = "Control")

keep <- rowSums(counts(dds)) > 2
dds <- dds[keep, ]
dds <- DESeq(dds,
  parallel = TRUE
)
print(resultsNames(dds))

res <- results(dds,
  name = "Group_Test_vs_Control",
  # contrast=c("Group","Test","Control")
)
summary(res)

# res2 为差异基因列表
res2 <- res %>%
  as_tibble(rownames = "ENSEMBL") %>%
  dplyr::left_join(geneMap, by = "ENSEMBL") %>%
  dplyr::select(SYMBOL, ENSEMBL, ENTREZID, everything()) %>%
  arrange(padj, desc(abs(log2FoldChange)))
write_delim(res2,
  file.path(
    outpath,
    paste0("DESeq2_DEG_", cohortname, ".txt")
  ),
  delim = "\t"
)

degFilter <- dplyr::filter(res2, abs(log2FoldChange) >= 1 & padj < 0.05) %>% 
  mutate(diff = if_else(log2FoldChange>0, 'up', 'down'))
write_delim(degFilter,
  file.path(
    outpath,
    paste0("DESeq2_DEG_filtered_", cohortname, ".txt")
  ),
  delim = "\t"
)


# DESeq2 提供了log2差异倍数(log2 fold change)收缩方法 lfcShrink
resLFC <- lfcShrink(dds,
  coef = "Group_Test_vs_Control",
  type = "apeglm",
  parallel = T
)
resLFC_df <- resLFC %>%
  as_tibble(rownames = "ENSEMBL") %>%
  dplyr::left_join(geneMap, by = "ENSEMBL") %>%
  dplyr::select(SYMBOL, ENSEMBL, ENTREZID, everything()) %>%
  arrange(padj)
write_delim(resLFC_df,
  file.path(
    outpath,
    paste0("DESeq2_DEG_lfcShrink_", cohortname, ".txt")
  ),
  delim = "\t"
)

# rlog vst,类似log2的转化，对counts进行转化。
# rldt 可用于下游的热图和聚类分析
ntd <- normTransform(dds) # this gives log2(n + 1)
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
rldt <- assay(rld) %>%
  as_tibble(rownames = "ENSEMBL") %>%
  dplyr::left_join(geneMap, by = "ENSEMBL") %>%
  dplyr::select(SYMBOL, ENSEMBL, ENTREZID, everything())

# 标准化的counts数据，normalized=F时，为原始counts数
# 计算量化因子, DESeq函数包装了四步
nReadCount <- counts(dds, normalized = TRUE) %>%
  as_tibble(rownames = "ENSEMBL") %>%
  dplyr::left_join(geneMap, by = "ENSEMBL") %>%
  dplyr::select(SYMBOL, ENSEMBL, ENTREZID, everything())


# basic plot by DESeq2 ----------------------------------------------------

plotMA(res, ylim = c(-2, 2))
plotMA(resLFC,
  ylim = c(-2, 2),
  alpha = 0.05
)

# 绘制某一个基因的counts
plotCounts(dds, gene = which.min(res$padj), intgroup = "Group")


# 对于绘图，其呈现的依旧是counts数据，选择log/vst/rlog
# 热图
top20 <- head(res2, n = 20)$ENSEMBL
anno <- as.data.frame(colData(dds)["Group"])

p_heat <- pheatmap::pheatmap(assay(ntd)[top20, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  clustering_method = "average",
  border = FALSE,
  annotation_col = anno
)

# PCA
p_pca <- DESeq2::plotPCA(vsd,
  intgroup = c("Group")
) +
  ggrepel::geom_label_repel(ggplot2::aes(label = name)) +
  ggplot2::theme_bw()

ggsave(file.path(
  outpath,
  paste0("PCA_", cohortname, ".png")
),
plot = p_pca, dpi = 800
)

pcaData <- plotPCA(vsd, intgroup = c("Group"), returnData = TRUE)
pc1 <- paste("PC1:", round(100 * attr(pcaData, "percentVar"))[1], "%")
pc2 <- paste("PC2:", round(100 * attr(pcaData, "percentVar"))[2], "%")

p_pca2 <- ggplot(pcaData, aes(PC1, PC2, label = name)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(values = c("red4", "green3", "orange", "blue2", "#C673FF2E")) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "transparent"),
    legend.title = element_blank(), legend.key = element_rect(fill = "transparent")
  ) +
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = "gray", size = 0.5) +
  geom_hline(yintercept = 0, color = "gray", size = 0.5) +
  ggrepel::geom_text_repel()


# 火山图
df <- res2 %>%
  mutate(logP = -log10(padj)) %>%
  drop_na(SYMBOL) %>%
  mutate(group = case_when(
    padj < 0.05 & log2FoldChange >= 1 ~ "up",
    padj < 0.05 & log2FoldChange <= -1 ~ "down",
    TRUE ~ "noSig"
  ))

df$label <- ""
up.genes <- head(df$SYMBOL[which(df$group == "up")], 10)
down.genes <- head(df$SYMBOL[which(df$group == "down")], 10)
top10.genes <- c(as.character(up.genes), as.character(down.genes))
df$label[match(top10.genes, df$SYMBOL)] <- top10.genes

p <- ggpubr::ggscatter(
  data = df,
  x = "log2FoldChange",
  y = "logP",
  color = "group",
  palette = c("#2f5688", "#BBBBBB", "#CC0000"),
  size = 1,
  label = df$label,
  font.label = 8,
  repel = T,
  xlab = "log2FoldChange",
  ylab = "-log10(adj)"
) + theme_bw() +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

ggsave(file.path(
  outpath,
  paste0("volcano_", cohortname, ".png")
),
plot = p, dpi = 800
)

# clusterprofiler ---------------------------------------------------------
# 差异基因
require(clusterProfiler)
require(org.Hs.eg.db)

gene_list <- degFilter$ENSEMBL

enrich.go <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL",
  pAdjustMethod = "BH",
  readable = TRUE,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
enrich.go@pvalueCutoff <- 1
enrich.go@qvalueCutoff <- 1

# enrich.go <- simplify(enrich.go)
pp <- clusterProfiler::dotplot(enrich.go,
  font.size = 8,
  color = "p.adjust",
  showCategory = 20,
  split = 'ONTOLOGY'
) + facet_grid(ONTOLOGY~., scale="free") +
  scale_y_discrete(labels = function(x) {
    str_wrap(x, width = 40)
  })

write_delim(enrich.go@result,
  file = file.path(
    outpath,
    paste0("go_result_", cohortname, ".txt")
  ),
  delim = "\t"
)
ggsave(file.path(
  outpath,
  paste0("go_dotplot_", cohortname, ".png")
),
plot = pp, dpi = 800
)

# KEGG
gene.Map <- clusterProfiler::bitr(gene_list,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = "org.Hs.eg.db"
)
kegg <- enrichKEGG(
  gene = genelist,
  keyType = "kegg",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
kegg@pvalueCutoff <- 1
kegg@qvalueCutoff <- 1

kegg_p <- clusterProfiler::dotplot(kegg,
  font.size = 8,
  showCategory = 20
) +
  scale_y_discrete(labels = function(x) {
    str_wrap(x, width = 40)
  })

kegg_p2 <- barplot(kegg,
  font.size = 8,
  showCategory = 20
) +
  scale_y_discrete(labels = function(x) {
    str_wrap(x, width = 40)
  })

ggsave(file.path(
  outpath,
  paste0("kegg_bapplot_", cohortname, ".png")
),
plot = kegg_p2, dpi = 800
)
ggsave(file.path(
  outpath,
  paste0("kegg_dotplot_", cohortname, ".png")
),
plot = kegg_p, dpi = 800
)
write_delim(kegg@result,
            file = file.path(
              outpath,
              paste0("kegg_result_", cohortname, ".txt")
            ),
            delim = "\t"
)

keggx <- setReadable(kegg)
p3 <- cnetplot(keggx, foldChange=geneList, 
               circular = TRUE, colorEdge = TRUE)

# gsea
gmt <- read.gmt('')
gene_sorted <- sort(gene_list)

gsea <- GSEA(
  geneList = gene_sorted,
  TERM2GENE = gmt
)

em2 <- setReadable(gsea, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")
em2.df <- as.data.frame(em2)
em2.df <- em2.df %>%
  dplyr::select(ID,NES,setSize,p.adjust)
GSEA_result2 <- merge(lipid_zz,em2.df, by.x='ont', by.y = 'ID')
id <- 1
anno <- em2[id, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

enrichplot::gseaplot2(gsea, 1) 
annotate("text", 0.8, 0.85,
         label = lab, hjust=0, vjust=0)

gse_go <- gseGO(
  geneList = gene_list,
  ont = 'ALL',#one of “BP”, “MF”, “CC” or “ALL”
  keyType = 'ENSEMBL',
  nPerm = 5000,
  minGSSize = 3, 
  maxGSSize = 800, 
  pvalueCutoff = 0.05, 
  verbose = TRUE, 
  OrgDb = organism, 
  pAdjustMethod = "none"
)



gseaplot(gse, by = "all", title = gse$Description[1], 
         geneSetID = 1
         )



# plot by FPKM ------------------------------------------------------------

fpkm <- read_delim(countfile, delim = "\t") %>%
  dplyr::select(`gene_id`, `Exonic.gene.sizes`, ends_with("FPKM")) %>%
  dplyr::select(`gene_id`, matches(paste(rownames(sampleGroup), collapse = "|"))) %>%
  column_to_rownames("gene_id")
colnames(fpkm) <- sapply(colnames(fpkm), function(x) {
  str_split(x, "_FPKM")[[1]][1]
})

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(fpkm,2,fpkmToTpm)

top20 <- head(res2, n = 30)$ENSEMBL
# top20 <- degFilter$ENSEMBL #按照差异基因绘制热图
anno <- as.data.frame(colData(dds)["Group"])
fpkm_top20 <- fpkm[top20, ] %>%
  rownames_to_column() %>%
  left_join(geneMap, by = c("rowname" = "ENSEMBL")) %>%
  dplyr::select(-c(rowname, ENTREZID)) %>%
  drop_na(SYMBOL) %>%
  column_to_rownames("SYMBOL")

p_heat <- pheatmap::pheatmap(log10(fpkm_top20 + 1),
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  clustering_method = "average",
  border = FALSE,
  annotation_col = anno
)
ggsave(file.path(
  outpath,
  paste0("heatmap_FPKM_", cohortname, ".png")
),
plot = p_heat, dpi = 800
)

# PCA
df3 <- fpkm %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame()
df3 <- log10(df3 + 1)
res.pca <- prcomp(df3, scale = FALSE)
df_out <- as.data.frame(res.pca$x)
percentage <- round(res.pca$sdev / sum(res.pca$sdev) * 100, 2)
percentage <- paste0(colnames(df_out), "(", paste(as.character(percentage), "%", ")", sep = ""))
pc1 <- percentage[1]
pc2 <- percentage[2]

df4 <- merge(df_out, anno, by = 0)

p_pca3 <- ggplot(df4, aes(PC1, PC2, label = `Row.names`)) +
  geom_point(aes(color = Group)) +
  scale_color_manual(values = c("red4", "green3", "orange", "blue2", "#C673FF2E")) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "transparent"),
    legend.title = element_blank(), legend.key = element_rect(fill = "transparent")
  ) +
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = "gray", size = 0.5) +
  geom_hline(yintercept = 0, color = "gray", size = 0.5) +
  ggrepel::geom_text_repel()

ggsave(file.path(
  outpath,
  paste0("PCA_FPKM_", cohortname, ".png")
),
plot = p_pca3, dpi = 800
)

