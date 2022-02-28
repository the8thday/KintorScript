#

library(GWENA)
library(magrittr) # Not mandatory, we use the pipe `%>%` to ease readability.

threads_to_use <- 2


# prepare data
expr <- '/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/re_analysis/group_12_paired/DESeq2_rld_group12.txt'

expr <- read_delim(expr, delim = '\t') %>%
  distinct() %>%
  drop_na() %>%
  dplyr::select(-c(ENTREZID,ENSEMBL)) %>%
  group_by(SYMBOL) %>% top_n(1, abs(`12506I190007`)) %>% ungroup() %>% tibble() %>%
  column_to_rownames('SYMBOL') %>% t() %>%
  as.data.frame()

metafile <- "/Users/congliu/OneDrive/kintor/Daily_Work/follicle_RNA/sample25.xlsx"
sampleGroup <- readxl::read_excel(metafile) %>%
  filter(class != 'normal')

is_data_expr(expr)

# Gene filtering
expr_filtered <- filter_low_var(expr, pct = 0.7, type = "median")

ncol(expr_filtered)

# Network building
net <- build_net(expr_filtered,
                 fit_cut_off = 0.9,
                 cor_func = "spearman",
                 n_threads = threads_to_use)
#Power selected :
net$metadata$power

# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

modules <- detect_modules(expr_filtered,
                          net$network,
                          detailled_result = TRUE,
                          merge_threshold = 0.25)

# Number of modules before merging :
length(unique(modules$modules_premerge))
# Number of modules after merging:
length(unique(modules$modules))

layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge,
  modules_merged = modules$modules)


# Functional enrichment
enrichment <- bio_enrich(modules$modules)
plot_enrichment(enrichment)

# Phenotypic association
# With data.frame/matrix
phenotype_association <- associate_phenotype(
  modules$modules_eigengenes,
  sampleGroup %>% dplyr::select(Condition, Age, Slide))
plot_modules_phenotype(phenotype_association)


# Graph visualization and topological analysis, get hub gene

module_example <- modules$modules$`2` # 依据modules after merging选择
# graph <- build_graph_from_sq_mat(net$network[module_example, module_example]) # 挺吃资源,不要运行

layout_mod_2 <- plot_module(graph, upper_weight_th = 0.999995,
                            vertex.label.cex = 0,
                            node_scaling_max = 7,
                            legend_cex = 1)
# detect the sub clusters inside a module to find genes working toward the same function through enrichment
net_mod_2 <- net$network[modules$modules$`2`, modules$modules$`2`]
sub_clusters <- get_sub_clusters(net_mod_2)

layout_mod_2_sub_clust <- plot_module(graph, upper_weight_th = 0.999995,
                                      groups = sub_clusters,
                                      vertex.label.cex = 0,
                                      node_scaling_max = 7,
                                      legend_cex = 1)

# hub gene
get_hub_high_co(network = net$network, # Highest connectivity
                modules = modules$modules,
                top_n = 5
                )

get_hub_degree(network = net, # Superior degree
               modules = modules
               )

get_hub_kleinberg(network = net, # Kleinberg’s score
                  modules = modules,
                  top_n = 5)







