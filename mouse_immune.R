# 小鼠组织免疫浸润分析


library(tidyverse)
source('/Users/congliu/OneDrive/kintor/Script/RNA_analysis/RNA_R_script/CIBERSORTx.R')



# fpkm as input -----------------------------------------------------------

fpkm <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/all.fpkm_anno.xls',
                   delim = '\t'
                   )
df <- fpkm %>% select(gene_id, ends_with('FPKM'))

# write_delim(df, file = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/counts.txt',
#             delim = '\t'
#             )
geneMap <- clusterProfiler::bitr(df$gene_id,
                                 fromType = "ENSEMBL",
                                 toType = c("SYMBOL", "ENTREZID"),
                                 OrgDb = "org.Mm.eg.db"
                                 )

foo <- df %>% left_join(geneMap, by = c('gene_id'='ENSEMBL')) %>%
  dplyr::select(-c(gene_id, ENTREZID)) %>%
  drop_na() %>% distinct() %>%
  group_by(SYMBOL) %>% top_n(1, abs(`G4-2_FPKM`)) %>% ungroup() %>%
  filter(!SYMBOL %in% c('Fam90a1b', 'Gm21470', 'Vmn2r122', 'Vmn2r29')) %>%
  dplyr::select(SYMBOL, everything())
# write_delim(foo, file = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/exp.file_CIBERSORT.txt',
#             delim = '\t')


# 准备CIBERSORT的输入文件
# sig_matrix <- '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/waip_sig_mice.txt'
sig_matrix <- ''
mixture_file <- '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/exp.file_CIBERSORT.txt'
sampleGroup <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/sampleInfor.xlsx')

res_cibersort = CIBERSORT(sig_matrix = sig_matrix,
                          mixture_file = mixture_file,
                          perm = 1000, QN = FALSE)

head(res_cibersort)
write_delim(data.frame(res_cibersort) %>% rownames_to_column(),
            file = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/CIBERSORT.txt',
            delim = '\t')

df_plot <- as.data.frame(res_cibersort) %>%
  rownames_to_column() %>%
  left_join(sampleGroup, by = c('rowname'='sample'))

df_plot1 <- df_plot %>% pivot_longer(cols = 2:26,
                                     names_to = 'CellType', values_to = 'Composition'
) %>% arrange(class)
df_plot1$rowname <- factor(df_plot1$rowname, levels = unique(df_plot1$rowname))

# boxplot
p1 <- ggpubr::ggboxplot(
  df_plot1,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition"
) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ),
  legend.position = 'none'
  )

ggsave(filename = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/CIBERSORT_boxpot.pdf',
       plot = p1,
       width = 10, height = 10
)



colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF",
            "#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C",
            "#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
            "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
            "#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="right",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_line(size=1),
        text = element_text(family="Times"),
        axis.text.y = element_text(size = 8,face='bold',color='black'),
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),
        axis.title = element_text(size=10,face="bold"),
        plot.title = element_text(hjust=0.5,size=10))
}


p2 <- ggplot(df_plot1, aes(rowname,Composition,fill = CellType)) +
  geom_bar(stat='identity') +
  # coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  my_theme()

ggsave(filename = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/CIBERSORT_bar2.pdf',
       plot = p2,
       width = 10, height = 10
       )




# counts ------------------------------------------------------------------

counts.raw <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/counts.txt',
                   delim = '\t'
                   ) %>% column_to_rownames('gene_id')
n <- nrow(counts.raw)
gene <- grep("ENS", rownames(counts.raw))
counts.raw <- counts.raw[gene, ]
load('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/receptor.ensemble.merge.RData')

counts.merge <- receptor_merge(counts.raw, gene.list=receptor.ensemble.merge)

p <- 0.75
counts.quatile <- quartile(counts.merge, p)

result.path <- '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/'
write.csv(counts.quatile,
          paste(result.path, "/ImmunCC.immuereceptor.quartileNorm.csv", sep=""),
          row.names=T, quote=F)

# CIBERSORT
mixture_file <- '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/ImmunCC.immuereceptor.quartileNorm.txt'
sig_matrix <- '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/SignatureMatrix.rnaseq.txt'
sampleGroup <- readxl::read_excel('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/sampleInfor2.xlsx')

res_cibersort = CIBERSORT(sig_matrix = sig_matrix,
                          mixture_file = mixture_file,
                          perm = 1000, QN = FALSE)

df_plot <- as.data.frame(res_cibersort) %>%
  rownames_to_column() %>%
  left_join(sampleGroup, by = c('rowname'='sample'))

df_plot1 <- df_plot %>% pivot_longer(cols = 2:11,
                                     names_to = 'CellType', values_to = 'Composition'
) %>% arrange(class)
df_plot1$rowname <- factor(df_plot1$rowname, levels = unique(df_plot1$rowname))


# mMCP-counter ------------------------------------------------------------

library("mMCPcounter")

fpkm <- read_delim('/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/all.fpkm_anno.xls',
                   delim = '\t'
)
df <- fpkm %>% select(gene_id, ends_with('FPKM'))

expressionData <- log2((df %>% column_to_rownames('gene_id')) + 1)

res_mMCP <- mMCPcounter.estimate(expressionData, features = "ENSEMBL.ID")

head(res_mMCP)
write_delim(data.frame(res_mMCP) %>% rownames_to_column(),
            file = '/Users/congliu/OneDrive/kintor/Daily_Work/小鼠组织RNA/mMCP_res.txt',
            delim = '\t')

df_plot <- as.data.frame(res_mMCP) %>% t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  left_join(sampleGroup, by = c('rowname'='sample'))

df_plot1 <- df_plot %>% pivot_longer(cols = 2:17,
                                     names_to = 'CellType', values_to = 'Composition'
) %>% arrange(class)
df_plot1$rowname <- factor(df_plot1$rowname, levels = unique(df_plot1$rowname))





