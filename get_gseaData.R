# 提取clusterProfiler中GSEA的结果

# gseaScores <- getFromNamespace("gseaScores", "DOSE")
gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES,
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList

  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]

  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df_l <- gseaScores(geneList, geneSet, exponent, fortify=FALSE)
  df <- df_l$runningES
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList

  df$Description <- object@result[geneSetID, "Description"]
  res <- list(ES = ES, runningES = df)
  return(res)
}

head(gsdata <- gsInfo(gsea,1)$runningES)
ES <- gsInfo(gsea,1)$ES

setGene <- subset(gsdata, position != 0)
leading_gene <- (gsea@result %>% pull(core_enrichment) %>%
                   str_split(pattern = '/'))[[1]]
my_gene <- setGene %>% filter(gene %in% leading_gene)
# 提取顶点基因
topgene <- tail(my_gene ,1)

# plot
p <- ggplot(gsdata, aes_(x = ~x, y = ~runningScore)) + xlab(NULL) +
  theme_classic(11) +
  theme(panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(expand=c(0,0)) +
  geom_line(aes_(color= ~Description),
            size=1) +
  theme(legend.position = c(.8, .8), legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))


p2 <- p + geom_point(data = topgene,aes(x=x,y=runningScore),
                      color='red',size=3) +
  geom_segment(data = topgene,aes(xend=x,yend=0),
               color='red',size=1) +
  geom_text_repel(data = my_gene,aes(x=x,y=runningScore,label=gene),
                     color='black') +
  geom_segment(data = setGene, aes(x = x, y = runningScore, color = geneList, xend=x,yend=0)) +
  scale_color_gradient2(mid = "green")


pp <- ggplot(data = setGene,
       aes(x = x, y = runningScore, color = geneList)) +
  geom_segment(aes(xend=x,yend=0)) +
  scale_color_gradient2(mid = "darkolivegreen4") +
  geom_line(size=1) +
  theme_classic() +
  ylab('Runing Enrichment score (ES)') +
  xlab('Position in the Ranked List of Genes') +
  theme(axis.title = element_text(size = 18,face = 'bold'),
        legend.position = 'none',
        axis.text = element_text(size = 16,face = 'bold'),
        axis.line = element_line(size = 1),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(size = 18,face = 'bold')
  ) +
  geom_hline(yintercept = 0,size=0.5, color = 'gray27') +
  ggtitle('GSEA Enrichment') +
  geom_point(data = topgene,aes(x=x,y=runningScore),
             color='firebrick3',size=2) +
  geom_segment(data = topgene,aes(xend=x,yend=0),
               color='firebrick3',size=1) +
  ggrepel::geom_text_repel(data = my_gene,
                           mapping = aes(x = x, y = runningScore, label=gene),
                  color='black',
                  force = 1,
                  max.overlaps = 10,
                  arrow = arrow(length = unit(0.02, "npc"),
                                ends = "last", type = "open"),
                  segment.color= 'black') +
  annotate('text', x = 10000, y = 0.75,
           label = gsea@result[1, "ID"]
           )
pp



