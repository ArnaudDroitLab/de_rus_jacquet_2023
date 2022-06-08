# Données téléchargées le 29/05/2022
# GSEA : https://www.gsea-msigdb.org/gsea/msigdb/cards/ST_ERK1_ERK2_MAPK_PATHWAY  (32 gènes)
# Biocarta : https://maayanlab.cloud/Harmonizome/gene_set/erk1%24slash%24erk2+mapk+signaling+pathway/Biocarta+Pathways  (20 gènes)
# Kegg : https://www.genome.jp/dbget-bin/www_bget?pathway:map04010 (229 gènes)
# QuickGO : https://www.ebi.ac.uk/QuickGO/annotations?goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:0070371&taxonId=9606&taxonUsage=descendants (56 gènes)

library(tidyverse)
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)

setwd("/mnt/tn02_bioinfo/JACA/IPS/post_processing/ERK_pathway_db")

quickGO <- read_tsv("QuickGO-annotations-1654105854851-20220601.tsv") %>%
  arrange(SYMBOL) %>%
  mutate_all(.funs = toupper) %>%
  pull(SYMBOL) %>%
  unique %>%
  mapIds(x = org.Hs.eg.db, keytype = "SYMBOL", column = "ENSEMBL") %>%
  as.data.frame %>%
  rownames_to_column("SYMBOL") %>%
  `colnames<-` (c("SYMBOL", "ENSEMBL"))

GSEA <- read_csv("GSEA.csv") %>%
  arrange(MAPPED_SYMBOLS) %>%
  mutate_all(.funs = toupper) %>%
  pull(MAPPED_SYMBOLS) %>%
  unique %>%
  mapIds(x = org.Hs.eg.db, keytype = "SYMBOL", column = "ENSEMBL") %>%
  as.data.frame %>%
  rownames_to_column("SYMBOL") %>%
  `colnames<-` (c("SYMBOL", "ENSEMBL"))

bioCarta <- read_csv("bioCarta.csv") %>%
  arrange(Symbol) %>%
  mutate_all(.funs = toupper) %>%
  pull(Symbol) %>%
  unique %>%
  mapIds(x = org.Hs.eg.db, keytype = "SYMBOL", column = "ENSEMBL") %>%
  as.data.frame %>%
  rownames_to_column("SYMBOL") %>%
  `colnames<-` (c("SYMBOL", "ENSEMBL"))
 
kegg <- read_csv("kegg.csv") %>%
  arrange(symbol) %>%
  mutate_all(.funs = toupper) %>%
  pull(symbol) %>%
  unique %>%
  mapIds(x = org.Hs.eg.db, keytype = "SYMBOL", column = "ENSEMBL") %>%
  as.data.frame %>%
  rownames_to_column("SYMBOL") %>%
  `colnames<-` (c("SYMBOL", "ENSEMBL"))

dim(quickGO)
dim(GSEA)
dim(bioCarta)
dim(kegg)

all_symbols_list <- list("quickGO" = quickGO,
                         "GSEA" = GSEA,
                         "bioCarta" = bioCarta,
                         "kegg" = kegg) %>%
  write.xlsx("ERK_pathway_genes_DB_to_filter.xlsx")

all_ensembl <- c(quickGO$ENSEMBL, GSEA$ENSEMBL, bioCarta$ENSEMBL, kegg$ENSEMBL) %>%
  table %>%
  as.data.frame %>%
  arrange(desc(Freq))
table(all_ensembl$Freq)


# Après vérification, à la main, des symbols dont les ensembl n'ont pas été identifiés
xlsx_quickGO <- read.xlsx("ERK_pathway_genes_DB.xlsx", sheet = 1)
xlsx_GSEA <- read.xlsx("ERK_pathway_genes_DB.xlsx", sheet = 2)
xlsx_bioCarta <- read.xlsx("ERK_pathway_genes_DB.xlsx", sheet = 3)
xlsx_kegg <- read.xlsx("ERK_pathway_genes_DB.xlsx", sheet = 4)

dim(quickGO)
dim(GSEA)
dim(bioCarta)
dim(kegg)

dim(xlsx_quickGO)
dim(xlsx_GSEA)
dim(xlsx_bioCarta)
dim(xlsx_kegg)

all_ensembl <- c(xlsx_quickGO$ENSEMBL, xlsx_GSEA$ENSEMBL, xlsx_bioCarta$ENSEMBL, xlsx_kegg$ENSEMBL) %>%
  table %>%
  as.data.frame %>%
  arrange(desc(Freq))
table(all_ensembl$Freq)
