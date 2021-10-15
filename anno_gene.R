# annotate single gene

library(biomaRt)
library(ensembldb)
library(mygene)


# first explore ensemble --------------------------------------------------

# look at top 10 databases      
head(biomaRt::listMarts(), 10)
listEnsembl() # list biomart
listEnsemblArchives()

# 主要是两个参数
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      # mirror = 'asia'
                      )

head(listFilters(ensembl)) # list what to find, what you have
listAttributes(ensembl) # list what you need
searchAttributes('hpa', mart = ensembl)

entrez=c("673","837")

getBM(
  attributes = c('ensembl_gene_id', 'go_id', 'hgnc_symbol',
                 'description', 'phenotype_description',
                 'kegg_enzyme', 'name_1006','hpa_accession','mim_gene_description'
                 ),
  filters = c('entrezgene_id'),
  values = entrez,
  mart = ensembl
)



# next explore mygene -----------------------------------------------------

genelist <- c('AR', 'MYC')

getGenes(
  geneids = entrez,
  species = 'human',
  fields = 'all'
)

queryMany(genelist, 
          scopes="symbol", 
          fields="entrezgene", 
          species="human")








