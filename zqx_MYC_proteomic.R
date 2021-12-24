## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2021-8-17
##
## Copyright (c) cliu, 2021
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: 为zqx所提供的蛋白质谱数据，求差异表达，并进行hub gene分析，网络图拓扑结构分析.
##         于12月份再次进行分析
##
## ---------------------------

library(tidyverse)


df_profiling_3nM <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/qianx/12-17/3nM_profiling_with significant.xlsx')
df_profiling_1.5nM <- readxl::read_excel(
  '/Users/congliu/OneDrive/kintor/qianx/12-17/1.5nM_profiling_with significant.xlsx')

# 参数
df <- df_profiling_3nM
path = '/Users/congliu/OneDrive/kintor/qianx/12-17'
cohert = '3nM'

df$fdr <- p.adjust(df$ttest, method = 'BH')

df_diff <- df %>% filter(signficant %in% c('up', 'down')) %>%
  arrange(fdr, desc(abs(log2ratio)))
dim(df_diff)
write_delim(df_diff,
            file.path(path, glue::glue('{cohert}_volcano_genelist.txt')),
            delim = '\t'
            )

df_p <- df %>%
  mutate(baseMean = log2(sqrt(mean_ctrl * mean_3nM)))

# 绘制火山图
top10_up <- df_p %>% filter(signficant == 'up') %>%
  arrange(fdr, desc(abs(log2ratio))) %>%
  slice_head(n=10)
top10_down <- df_diff %>% filter(signficant == 'down') %>%
  arrange(fdr, desc(abs(log2ratio))) %>%
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
  labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)', color = '', title = glue::glue('{cohert} vs DMSO')) +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = log2ratio, y = -log10pvalue,
                                                                       label = GeneSymbol),
                           size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black',
                           show.legend = FALSE)

ggsave(plot = p, filename = file.path(path, glue::glue('{cohert}_volcano.pdf')),
       width = 10, height = 10)

# MAplot
p2 <- ggplot(df_p, aes(x = baseMean, y = log2ratio, color = signficant)) +
  geom_point(size = 2) +
  scale_colour_manual(values  = c("#B31B21", "#1465AC", "darkgray"), limits = c('up', 'down', 'no')) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'), legend.position = 'top') +
  labs(x = '\nLog2 Mean Expression', y = 'Log2 Fold Change\n', color = '') +
  geom_hline(yintercept = c(0, -log2(2), log2(2)), linetype = c(1, 2, 2),
             color = c("black", "black", "black")) +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = baseMean, y = log2ratio,
                                                                       label = GeneSymbol),
                           size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)

ggsave(plot = p2, filename = file.path(path, glue::glue('{cohert}_MAplot.pdf')),
       width = 10, height = 10)



# string
library(STRINGdb)
library(igraph)
library(ggraph)

gene.df <- df_diff %>% select(GeneSymbol) %>%
  as.data.frame()

string_db <- STRINGdb$new(species=9606,
                          version='11.5',
                          score_threshold=400
                          )
STRINGdb$methods()

data_mapped <- string_db$map(my_data_frame = gene.df,
                             my_data_frame_id_col_names = "GeneSymbol",
                             takeFirst = TRUE,
                             removeUnmappedRows = TRUE)
print(dim(data_mapped))
string_db$plot_network(data_mapped$STRING_id)
data_links <- data_mapped$STRING_id %>% string_db$get_interactions()
head(data_links)

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
head(net)

igraph::V(net)$deg <- igraph::degree(net)
igraph::V(net)$size <- igraph::degree(net)/5
igraph::E(net)$width <- igraph::E(net)$weight/10
# 使用ggraph绘图
ggraph(net,layout = "stress")+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = T)+
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

ggsave(plot = p2, filename = file.path(path, glue::glue('{cohert}_PPI.pdf')),
       width = 10, height = 10)


# eulerr ------------------------------------------------------------------

library(eulerr)
library(ggvenn)

diff_3nm <- read_delim(file.path(path, '3nM_volcano_genelist.txt'), delim = '\t') %>%
  pull(GeneSymbol)
diff_1.5nm <- read_delim(file.path(path, '1.5nM_volcano_genelist.txt'), delim = '\t') %>%
  pull(GeneSymbol)


a <- list(nM3 = diff_3nm,
          nM1.5 = diff_1.5nm)

ggvenn(
  a,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggvenn(a, show_elements = TRUE, label_sep = "\n")




