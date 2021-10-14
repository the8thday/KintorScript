library(tidyverse)


df <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/1.5nM_profiling_with significant.xlsx')

df %>% filter(ttest <= 0.05) %>%
  filter(abs(log2ratio)>=1)

path = '/Users/congliu/OneDrive/kintor/qianx/'
cohert = '1.5nM'


df_p <- df %>% filter(signficant %in% c('up', 'down')) %>%
  arrange(desc(abs(log2ratio)), ttest)

i_gene <- Hmisc::Cs(
  GSPT2,
  GSPT1,
  CK1,
  Ikaros,
  Aiolos
) %>% str_to_upper()


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
  geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
  # xlim(-5, 5) +
  # ylim(0, 6) +
  labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)', color = '', title = '1.5nM vs DMSO') +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = log2ratio, y = -log10pvalue,
                                                                       label = GeneSymbol),
                           size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black',
                           show.legend = FALSE)

ggsave(plot = p, filename = file.path(path, glue::glue('{cohert}_volcano.pdf')),
       width = 10, height = 10, dpi = 500)









