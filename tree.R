
library(tidyverse)
library(ggtree)
library(ggtreeExtra)


tree = read.tree('d:/Brazil_seq/h348_1611/CovMutation/cov100/cov100_aln.treefile')
as_tibble(tree)

# group 信息视情况而定
groupInfo <- read.delim('~/Downloads/sample_group.txt', sep = '\t', header = T,
                        row.names = 1
)

tree_g <- groupOTU(tree, groupInfo)


p <- ggtree(tree, layout="rectangular", 
            size=0.5, col="deepskyblue4") + 
  geom_tiplab(size=1, color="gray15", align=T, 
              linetype=3, linesize=0.5, hjust=-0.02) +
  geom_tippoint(size=1.5, color="deepskyblue4") + 
  # geom_text2(aes(subset=!isTip, label=node), hjust=-0.3, size=2, color="black") + 
  geom_nodepoint(color="orange", alpha=1/2, size=1) + 
  theme_tree2()

p2 <- ggtree(tree, layout="circular", ladderize=FALSE, size=0.8, 
                branch.length="none") +
  geom_tiplab2(hjust=-0.3) +
  geom_tippoint(size=3)+
  theme(legend.title=element_text(face="bold"), legend.position="bottom", 
        legend.box="horizontal", legend.text=element_text(size=rel(0.8)))

