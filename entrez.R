# 批量下载文献


# Entrez Utilities 为其提供的api接口

library(rentrez)

# show all avilible database
entrez_dbs()

# what the database is
entrez_db_summary('pubmed')


# 看某个数据库的检索关键词
# search field is very important!
entrez_db_searchable("pubmed")


res <- entrez_search(
  db = 'pubmed',
  term = 'TP53 & prostate cancer AND 2015:2021[PDAT]'
)

summary(res)
res$ids












