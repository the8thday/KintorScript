# use rentrez to download ncbi info

library(rentrez)

# 所有数据库
entrez_dbs()

# 看基本信息
entrez_db_summary("pubmed")


# 看某个数据库的检索关键词
entrez_db_searchable("sra")


r_search <- entrez_search(db="pubmed", term="single cell[All Fields]")
r_search

# 加日期限定
entrez_search(db="sra",
              term="Homo sapiens[ORGN] AND 2019:2020[PDAT]",
              retmax=0)



# 找到与某个基因有关联的其他基因

gene_search <- entrez_search(db="gene", term="(ZNF365[GENE]) AND (Homo sapiens[ORGN])")
gene_search
interactions_from_gene <- function(gene_id){
  xmlrec <- entrez_fetch(db="gene", id=gene_id, rettype="xml", parsed=TRUE)
  XML::xpathSApply(xmlrec,
                   "//Gene-commentary[Gene-commentary_heading[./text()='Interactions']]//Other-source[Other-source_src/Dbtag/Dbtag_db[./text()='GeneID']]//Other-source_anchor",
                   XML::xmlValue)
}
interactions_from_gene(gene_search$ids)

# 根据核酸记录ID，提取相关文献

library(rentrez)
library(XML)
library(dplyr)

# function to add column to df if not already included
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- NA
  data
}

# function to collate all publications associated with sequences
get_pub_info <- function(i){
  fetch2 <- entrez_fetch(db = "nucleotide", id = i,
                         rettype = "gbc", retmode="xml", parsed = TRUE)
  xml_list2 <- xmlToList(fetch2)
  ref_list <- xml_list2$INSDSeq$INSDSeq_references
  # extract publication fields info
  authors <- unlist(ref_list$INSDReference$INSDReference_authors) %>% paste(collapse = "; ")
  title <- ref_list$INSDReference$INSDReference_title
  journal <- ref_list$INSDReference$INSDReference_journal
  year <-gsub(".*\\((.*)\\).*", "\\1", journal)
  pm_id <- ref_list$INSDReference$INSDReference_pubmed
  remark <- ref_list$INSDReference$INSDReference_remark
  # create data frame of information
  pub.data <- data.frame(i, authors, journal, year)
  if(is.null(title)==FALSE) pub.data$title <- title
  if(is.null(pm_id)==FALSE) pub.data$pubmed_id <- pm_id
  if(is.null(remark)==FALSE) pub.data$remark <- remark
  pub.data <- fncols(pub.data, c("title", "pubmed_id", "remark"))
}

sequence_list <- c("AB687721.2", "AB600942.1", "AJ880277.1")
list_of_dfs <- lapply(sequence_list, get_pub_info) # run function on list of sequences
df_combine <- bind_rows(list_of_dfs)
colnames(df_combine)[1] <- "NCBI_idv"
df_combine <- tidyr::separate(df_combine, remark, c("text", "doi"), sep = "DOI:") # extract doi
df_combine <- tidyr::separate(df_combine, doi, c("doi", "text2"), sep = ";")
df_combine$remark <- paste(df_combine$text,df_combine$text2)
df_combine$text <- NULL
df_combine$text2 <- NULL

head(df_combine[, c("pubmed_id", "doi", "journal")])




