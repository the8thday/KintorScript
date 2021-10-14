library(tidyverse)


df <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/1.5nM_profiling_with significant.xlsx')

df_profiling_3nM <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/MYC_profiling.xlsx') %>%
  mutate(log2ratio = log2(ratio_3nM_ctrl),
         log10pvalue = log10(ttest_3nM_DMSO),
         signficant = case_when(
           log2ratio>=1 & ttest_3nM_DMSO<=0.05 ~ 'up',
           log2ratio<=-1 & ttest_3nM_DMSO<=0.05 ~ 'down',
           TRUE ~ 'no'
         )
         ) %>%
  rename('ttest'='ttest_3nM_DMSO')
df_profiling_3nM$fdr <- p.adjust(df_profiling_3nM$ttest, method = 'BH')

fdr <- p.adjust(df$ttest, method = 'BH')
df$fdr <- fdr

path = '/Users/congliu/OneDrive/kintor/qianx/'
cohert = '3nM'

df <- df_profiling_3nM
df_p <- df %>% filter(signficant %in% c('up', 'down')) %>%
  arrange(desc(abs(log2ratio)), ttest)
dim(df_p)
write_delim(df_p,
            file.path(path, glue::glue('{cohert}_volcano_gemelist.txt')),
            delim = '\t'
            )

i_gene <- Hmisc::Cs(
  GSPT2,
  GSPT1,
  CK1,
  Ikaros,
  Aiolos
) %>% str_to_upper()

df %>% filter(GeneSymbol %in% i_gene) # 并没有这几个基因的变化

# 绘制火山图
top10_up <- df_p %>% filter(signficant == 'up') %>%
  slice_head(n=10)
top10_down <- df_p %>% filter(signficant == 'down') %>%
  slice_head(n=10)

p <- ggplot(df, aes(x = log2ratio, y = -log10pvalue, color = signficant)) +
  geom_point(size = 2) +
  scale_colour_manual(values  = c("#B31B21", "#1465AC", "darkgray"), limits = c('up', 'down', 'no')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'),
        legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  # xlim(-5, 5) +
  # ylim(0, 6) +
  labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)', color = '', title = glue::glue('{cohert} vs DMSO')) +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = log2ratio, y = -log10pvalue,
                                                                       label = GeneSymbol),
                           size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black',
                           show.legend = FALSE)

ggsave(plot = p, filename = file.path(path, glue::glue('{cohert}_volcano.pdf')),
       width = 10, height = 10, dpi = 500)









