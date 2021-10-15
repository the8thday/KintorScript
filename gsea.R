#! /urs/bin/env Rscript

library(tidyverse)
# library(VariantAnnotation)
library(clusterProfiler)
# library(ChIPseeker)
library(org.Hs.eg.db)



setwd('/Node7Data/Data/g2/work/liucong/COV2/3.moabs/liuc_result/cov_pos_neg')

all_dmc <- read_delim('dmc_M2_cov_pos_2.G.bed_vs_neg_2.G.bed.txt', delim = '\t')
homer <- read_delim('all_dmc_covpos_covneg.xls', delim = '\t') %>% filter()#filter(abs(`Distance to TSS`) <= 2000)
all_dmc_key <- all_dmc %>% unite('key', `#chrom`:end)
homer_key <- homer %>% mutate(Start = Start - 1) %>% unite('key', Chr:End)
f_df <- homer_key %>% inner_join(all_dmc_key, by=c('key'='key')) %>% select(c(key, `Entrez ID`, `credibleDif_1-0`, `p_fet_1_v_0`, `Gene Name`, class))
write_delim(f_df, 'gsea_input_negnormal_005.txt', delim = '\t')


# for cpg gsea ------------------------------------------------------------

lipid_gmt <- '/Node7Data/Data/g2/work/liucong/COV2/3.moabs/liuc_result/cov_pos_neg/gmt_lipid/lipid_gene136_2k_cpg.gmt'
lipid_zz <- read.gmt(lipid_gmt)

# all_cpg <- read_delim('dmc_M2_cov_pos_2.G.bed_vs_neg_2.G.bed.txt', delim = '\t') %>% unite('key', `#chrom`:end) %>% 
#   arrange(desc(`credibleDif_1-0`)) %>% 
#   select(c(key, `credibleDif_1-0`, `p_fet_1_v_0`))
all_cpg <- read_delim('covpos_normal_for_gsea_005.txt', delim = '\t') %>% 
  arrange(desc(`credibleDif_1-0`))

df_list <- all_cpg$`credibleDif_1-0`
names(df_list) <- all_cpg$`#chrom_start_end`

gsea <- GSEA(df_list, TERM2GENE = lipid_zz, pvalueCutoff = 1)

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



# a intersect -------------------------------------------------------------


fisher0.05 <- all_dmc %>% filter(`p_fet_1_v_0` <= 0.05) %>% unite('key', `#chrom`:end)
credible_0.2_dmc <- read_delim('dmc_M2_cov_pos.G.bed_vs_normol_2.G.bed.bed', delim = '\t', skip = 1, col_names = F) %>% unite('key', X1:X3)
credible_0.2_dmc_2 <- all_dmc %>% filter(abs(`credibleDif_1-0`)>=0.2) %>% unite('key', `#chrom`:end)
credible_0.1_dmc_2 <- all_dmc %>% filter(abs(`credibleDif_1-0`)>=0.1) %>% unite('key', `#chrom`:end)
credible_0.05_dmc_2 <- all_dmc %>% filter(abs(`credibleDif_1-0`)>=0.05) %>% unite('key', `#chrom`:end)
credible_0.01_dmc_2 <- all_dmc %>% filter(abs(`credibleDif_1-0`)>=0.01) %>% unite('key', `#chrom`:end)

# intersect

fisher0.05_key <- fisher0.05$key
credible_0.2_dmc_2_key <- credible_0.2_dmc_2$key
credible_0.1_dmc_2_key <- credible_0.1_dmc_2$key
credible_0.05_dmc_2_key <- credible_0.05_dmc_2$key

length(intersect(credible_0.01_dmc_2$key, fisher0.05_key))

venn_list <- list(
  fisher0.05_key = fisher0.05$key,
  credible_0.2_dmc_2_key = credible_0.2_dmc_2$key,
  credible_0.1_dmc_2_key = credible_0.1_dmc_2$key,
  credible_0.05_dmc_2_key = credible_0.05_dmc_2$key,
  credible_0.01_dmc_2_key = credible_0.01_dmc_2$key
)

venn.plot <- venn.diagram(venn_list, filename = 'venn.tiff', 
                          fill = c('red', 'green', 'blue', 'yellow', 'pink'), alpha = 0.50,
                          # col = 'black', 
                          col = "transparent",
                          print.mode=c("raw"),
                          output=TRUE,
)







