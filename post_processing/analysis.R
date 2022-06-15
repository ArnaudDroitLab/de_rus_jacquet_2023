setwd("/mnt/tn02_bioinfo/JACA/IPS/post_processing/")
.libPaths(c("/home/cjbparlant/R/vls137", .libPaths()))
library(tidyverse)
library(rnaseq)
library(openxlsx)
library(factoextra)
library(ggrepel)
library(sva)
library("org.Hs.eg.db")
library(gridExtra)
library(ComplexHeatmap)

# Import
design <- read_csv("design.csv")

filenames <- rnaseq::get_filenames("../rnaseq_pipeline/results/kallisto/")
filenames <- filenames[which(names(filenames) %in% design$sample)]

anno <- "/is3/projects/PUBLIC/rnaseq_anno/org/Hs/Hs.Ensembl104.csv"

# txi <- import_kallisto(filenames, anno = anno)
# saveRDS(txi, "r_objects/txi.rds")
txi <- readRDS("r_objects/txi.rds")

raw_counts <- get_anno_df(txi, "raw_counts") %>%
  write.xlsx("livrables/raw_counts.xlsx")
tpm <- get_anno_df(txi, "tpm") %>%
  write.xlsx("livrables/tpm.xlsx")

# ------------------------------------ PCA ----------------------------------------------
# Produce PCA
pdf("livrables/pca.pdf")
res_pca <- produce_pca(txi)
dev.off()

res_pca_df <- res_pca$df %>%
  left_join(design)

# Plot PCA colored
p <- ggplot2::ggplot(res_pca_df, ggplot2::aes_string(x = "Dim1", y = "Dim2", color = "group", shape = "study")) +
     ggplot2::geom_point(size = 3) +
     ggrepel::geom_text_repel(ggplot2::aes(label = sample), color = "black", force = 10) +
     ggplot2::theme_bw() +
     ggplot2::xlab(paste0("Dim1 (", res_pca$pca$eig[1, 2] %>% round(2), "%)")) +
     ggplot2::ylab(paste0("Dim2 (", res_pca$pca$eig[2, 2] %>% round(2), "%)"))

pdf("livrables/pca_colored.pdf")
print(p)
dev.off()

# PCA Hierarchical clustering
res.hcpc <- FactoMineR::HCPC(res_pca$pca, graph = FALSE)
pdf("livrables/hierarchical_clustering.pdf")
fviz_dend(res.hcpc, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE,
           rect_border = "jco", labels_track_height = 0.8)
dev.off()

# ------------------------------------ PCA correction batch effect ----------------------------------------------
txi_combat <- txi
txi_combat$abundance <- ComBat_seq(txi_combat$abundance, batch=design$study, group=design$group)

pdf("livrables/pca_Combat_seq_colored.pdf")
res_pca_combat <- produce_pca(txi_combat)
dev.off()

res_pca_combat_df <- res_pca_combat$df %>%
  left_join(design)

p <- ggplot2::ggplot(res_pca_combat_df, ggplot2::aes_string(x = "Dim1", y = "Dim2", color = "group", shape = "study")) +
     ggplot2::geom_point(size = 3) +
     ggrepel::geom_text_repel(ggplot2::aes(label = sample), color = "black", force = 10) +
     ggplot2::theme_bw() +
     ggplot2::xlab(paste0("Dim1 (", res_pca_combat$pca$eig[1, 2] %>% round(2), "%)")) +
     ggplot2::ylab(paste0("Dim2 (", res_pca_combat$pca$eig[2, 2] %>% round(2), "%)"))

pdf("livrables/pca_Combat_seq_colored.pdf")
print(p)
dev.off()

res_hcpc_combat <- FactoMineR::HCPC(res_pca_combat$pca, graph = FALSE)
pdf("livrables/hierarchical_clustering_Combat_seq.pdf")
fviz_dend(res_hcpc_combat, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE,
           rect_border = "jco", labels_track_height = 0.8)
dev.off()

# ------------------------------------ DE ----------------------------------------------
dds <- deseq2_analysis(txi, design, ~ group)
saveRDS(dds, "r_objects/dds.rds")
dds <- readRDS("r_objects/dds.rds")

generate_DE_results <- function(group1, group2){
  res <- DESeq2::results(dds, contrast = c("group", group1, group2)) %>%
    as.data.frame %>%
    rownames_to_column("id") %>%
    left_join(txi$anno, by = "id") %>%
    dplyr::select(-id) %>%
    dplyr::select(ensembl_gene:transcript_type, everything()) %>%
    arrange(padj)
  write.xlsx(res, paste("livrables/de_", group1, "_vs_", group2, ".xlsx", sep = ""))

  # Volcano plots
  pdf(paste("livrables/volcano", group1, "vs", group2, "FC_1_5.pdf", sep = "_"))
  produce_volcano(res, fc_threshold = 1.5)$p +
    geom_text_repel(data = head(res, 10), mapping = aes(x = log2FoldChange, y = -log10(padj), label = symbol), color = "black")
  dev.off()

  pdf(paste("livrables/volcano", group1, "vs", group2, "FC_3.pdf", sep = "_"))
  produce_volcano(res, fc_threshold = 3)$p +
    geom_text_repel(data = head(res, 10), mapping = aes(x = log2FoldChange, y = -log10(padj), label = symbol), color = "black")
  dev.off()

  up_regulated <- res[res$log2FoldChange <= log2(1/3) & res$padj <= 0.05,] %>%
    na.omit() %>%
    write.xlsx(paste("livrables/up_regulated_", group1, "_vs_", group2, ".xlsx", sep = ""))

  down_regulated <- res[res$log2FoldChange >= log2(3) & res$padj <= 0.05,] %>%
    na.omit() %>%
    write.xlsx( paste("livrables/down_regulated_", group1, "_vs_", group2, ".xlsx", sep = ""))

  res
}

res_astro_IPSC <- generate_DE_results("Cortical_astrocytes", "Control_iPSCs")
res_astro_IPSCD_der <- generate_DE_results("Cortical_astrocytes", "Control_iPSC_derived_astrocytes")

# ------------------------------------ PCA normalisation VSD ----------------------------------------------
vsd <- vst(dds, blind=FALSE)

# Heatmap sample to sample distance
library("RColorBrewer")
library("pheatmap")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hm <- pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colors)

pdf("livrables/VSD_sample_to_sample_distance.pdf")
print(hm)
dev.off()

# PCA
pdf("livrables/VSD_PCA.pdf")
plotPCA(vsd, intgroup=c("study", "group"))
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("study", "group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("livrables/VSD_PCA_colored.pdf")
ggplot(pcaData, aes(PC1, PC2, color=study, shape=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()

# ------------------------------------ Facteurs de Transcription ----------------------------------------------
de_res_files <- dir(path = "livrables/livrables_20220330/DE", pattern = "de_", full.names = TRUE)
names(de_res_files) <- str_remove(de_res_files, "livrables/livrables_20220330/DE/de_") %>%
  str_remove(".xlsx")

de_res_list <- map(de_res_files, ~ read.xlsx(.x))

regulated_files <- dir(path = "livrables", pattern = "regulated", full.names = TRUE) # CF lines 109 to 115
# regulated_files <- regulated_files[!str_detect(regulated_files, "TF")]
names(regulated_files) <- str_remove(regulated_files, "livrables/") %>%
  str_remove(".xlsx")

regulated_list <- map(regulated_files, ~ read.xlsx(.x))

# DB 1 : The human Trancripition Factors db
TheHumanTrancriptionFactors_db <- read.csv("TF_db/TheHumanTrancriptionFactors_db.csv") %>%
  dplyr::select(Ensembl.ID)

# DB 2 : Human TFDB
HumanTFDB_db <- read_tsv("TF_db/HumanTFDB_db.txt") %>%
  dplyr::select(Ensembl)

# DB 3 : TF2DNA
pscan_files <- dir(path = "TF_db/TF2DNA_pscan_files", pattern = "pscan", full.names = TRUE, recursive = TRUE)

concat_all_pscan <- function(pscan_file, i, len){
  print(paste(i, "/", len))
  pscan_db <- str_remove(pscan_file, "TF_db/TF2DNA_pscan_files/Homo-sapiens_") %>%
    str_remove("/.*$")

  pscan_df <- read_tsv(pscan_file, show_col_types = FALSE) %>%
    mutate(database = pscan_db) %>%
    dplyr::select(tf_name, target_name, binding_score, p_value, database) %>%
    transform(tf_name = as.character(tf_name))
}

# TF2DNA_db <- map_dfr(pscan_files, ~ concat_all_pscan(.x, which(pscan_files == .x), length(pscan_files)))
# TF2DNA_db <- TF2DNA_db %>% write_csv("TF2DNA_db.csv")
TF2DNA_db <- read_csv("TF_db/TF2DNA_db.csv")

# Join all DB
join_TF_regulated_genes <- function(genes_df, TheHumanTrancriptionFactors_db, HumanTFDB_db, TF2DNA_db, name){
  print(name)
  genes_df$TheHumanTrancriptionFactors_db <- ""
  genes_df$TheHumanTrancriptionFactors_db[genes_df$ensembl_gene %in% TheHumanTrancriptionFactors_db$Ensembl.ID] <- "yes"
  genes_df$HumanTFDB_db <- ""
  genes_df$HumanTFDB_db[genes_df$ensembl_gene %in% HumanTFDB_db$Ensembl] <- "yes"
  genes_df$TF2DNA_db <- ""
  genes_df$TF2DNA_db[genes_df$symbol %in% unique(TF2DNA_db$tf_name)] <- "yes"
  write.xlsx(genes_df, paste("livrables/TF_", name, ".xlsx", sep = ""))
  genes_df
}

de_res_TF_list <- imap(de_res_list, ~ join_TF_regulated_genes(.x, TheHumanTrancriptionFactors_db, HumanTFDB_db, TF2DNA_db, .y))

# Add target columns
threshold <- quantile(TF2DNA_db$binding_score)[["75%"]]
TF2DNA_filtered <- TF2DNA_db %>% dplyr::filter(p_value < 0.0001, binding_score > threshold)

# filtrer sur les TF significatifs
get_TF_targets <- function(de_res_TF_df, de_name, regulated_list, TF2DNA_filtered){
  regulated_names <- names(regulated_list)[str_detect(names(regulated_list), de_name)]
  up <- regulated_list[[regulated_names[str_detect(regulated_names, "up")]]]
  down <- regulated_list[[regulated_names[str_detect(regulated_names, "down")]]]

  up_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% up$symbol) %>%
  dplyr::select(tf_name, target_name) %>%
  unique()

  down_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% down$symbol) %>%
  dplyr::select(tf_name, target_name) %>%
  unique()

  targeted_up_TF <- group_by(up_TF_targets, target_name) %>%
    summarize(targeted_by_up_regulated_TF  = paste(tf_name, collapse = ","))

  targeted_down_TF <- group_by(down_TF_targets, target_name) %>%
    summarize(targeted_by_down_regulated_TF = paste(tf_name, collapse = ","))

  de_res_TF_target_df <- left_join(de_res_TF_df, targeted_up_TF, by = c("symbol" = "target_name")) %>%
    left_join(targeted_down_TF, by = c("symbol" = "target_name"))

  all_de_res_TF_target_df <- list("de_res_TF_target" = de_res_TF_target_df,
                                  "targeted_up_TF" = targeted_up_TF,
                                  "targeted_down_TF" = targeted_down_TF)
  write.xlsx(all_de_res_TF_target_df, paste0("TF_target_", de_name, ".xlsx"))
}

TF_target_list <- imap(de_res_TF_list, ~ get_TF_targets(.x, .y, regulated_list, TF2DNA_filtered))


# ------------------------------------ Analyse ERK pathway ----------------------------------------------
# Get ERK pathway genes
ERK_genes <- read_csv("ERK_pathway_db/ERK_pathway_genes_allDB.csv") %>%
  dplyr::filter(ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$ensembl_gene) %>%
  dplyr::filter(SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$symbol) %>%
  dplyr::filter(ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$ensembl_gene) %>%
  dplyr::filter(SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$symbol)

stopifnot(ERK_genes$ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$ensembl_gene) # Check IDs
stopifnot(ERK_genes$SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$symbol) # Check IDs
stopifnot(ERK_genes$ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$ensembl_gene) # Check IDs
stopifnot(ERK_genes$SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$symbol) # Check IDs

# Get Figure 1 genes
fig1_genes <- read.xlsx("ERK_pathway_db/LRRK2_all datasets compiled_FINAL.xlsx", "List for TF_ERK analysis", colNames = FALSE) %>%
  unique %>%
  pull(X1) %>%
  mapIds(x = org.Hs.eg.db, keytype = "ENSEMBL", column = "SYMBOL") %>%
  as.data.frame %>%
  rownames_to_column("ENSEMBL") %>%
  `colnames<-` (c("ENSEMBL", "SYMBOL")) %>%
  arrange(SYMBOL)

stopifnot(fig1_genes$ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$ensembl_gene) # Check IDs
stopifnot(fig1_genes$SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSCs$symbol) # Check IDs
stopifnot(fig1_genes$ENSEMBL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$ensembl_gene) # Check IDs
stopifnot(fig1_genes$SYMBOL %in% de_res_TF_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$symbol) # Check IDs

get_TF_targets <- function(de_res_TF_df, de_name, regulated_list, TF2DNA_filtered, ERK_genes, fig1_genes){
  # up and down regulated genes analysis
  regulated_names <- names(regulated_list)[str_detect(names(regulated_list), de_name)] # get design name to extract regulated genes
  up <- regulated_list[[regulated_names[str_detect(regulated_names, "up")]]] # get up-regulated genes
  down <- regulated_list[[regulated_names[str_detect(regulated_names, "down")]]] # get down-regulated genes

  up_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% up$symbol) %>% # extract up-regulated TF
    dplyr::select(tf_name, target_name) %>%
    unique()
  targeted_up_TF <- group_by(up_TF_targets, target_name) %>% # extraction of all up-regulated TF targeting a gene
    summarize(targeted_by_up_regulated_TF  = paste(tf_name, collapse = ","))

  down_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% down$symbol) %>% # extract down-regulated TF
    dplyr::select(tf_name, target_name) %>%
    unique()
  targeted_down_TF <- group_by(down_TF_targets, target_name) %>% # extraction of all down-regulated TF targeting a gene
    summarize(targeted_by_down_regulated_TF = paste(tf_name, collapse = ","))

  # ERK pathway genes
  de_res_TF_df$is_ERK <- "" # Add column to identify genes regulated in ERK pathway
  de_res_TF_df$is_ERK[de_res_TF_df$ensembl_gene %in% ERK_genes$ENSEMBL] <- "yes"

  ERK_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% ERK_genes$SYMBOL) %>%
    dplyr::select(tf_name, target_name) %>%
    unique()
  targeted_ERK_TF <- group_by(ERK_TF_targets, target_name) %>%
    summarize(targeted_by_ERK_TF  = paste(tf_name, collapse = ","))

  de_res_TF_df$is_fig1 <- "" # Add column to identify regulated genes from Fig1
  de_res_TF_df$is_fig1[de_res_TF_df$ensembl_gene %in% fig1_genes$ENSEMBL] <- "yes"

  fig1_TF_targets <- dplyr::filter(TF2DNA_filtered, tf_name %in% fig1_genes$SYMBOL) %>%
    dplyr::select(tf_name, target_name) %>%
    unique()
  targeted_fig1_TF <- group_by(fig1_TF_targets, target_name) %>%
    summarize(targeted_by_fig1_TF  = paste(tf_name, collapse = ","))

  # Add new columns
  de_res_TF_target_df <- left_join(de_res_TF_df, targeted_up_TF, by = c("symbol" = "target_name")) %>% # add up regulated column
    left_join(targeted_down_TF, by = c("symbol" = "target_name")) %>% # add down-regulated column
    left_join(targeted_ERK_TF, by = c("symbol" = "target_name")) %>% # add ERK column
    left_join(targeted_fig1_TF, by = c("symbol" = "target_name")) # add fig1 column

  # Write xlsx output
  all_de_res_TF_target_df <- list("de_res_TF_target" = de_res_TF_target_df,
                                  "up_TF_targets" = up_TF_targets, # add sheet with genes targeted by up-regulated TF
                                  "down_TF_targets" = down_TF_targets, # add sheet with genes targeted by up-regulated TF
                                  "ERK_TF_targets" = ERK_TF_targets, # add sheet with genes targeted by up-regulated TF
                                  "fig1_TF_targets" = fig1_TF_targets) # add sheet with genes targeted by up-regulated TF
  write.xlsx(all_de_res_TF_target_df, paste0("TF_target_ERK_", de_name, ".xlsx"))

  de_res_TF_target_df
}

TF_target_list <- imap(de_res_TF_list, ~ get_TF_targets(.x, .y, regulated_list, TF2DNA_filtered, ERK_genes, fig1_genes))

# ------------------------------------ Visualisation des résults ----------------------------------------------
## Histogrammes
# Histogramme du nombres de gènes de la fig1 ciblés par les FT ERK comparés au nombres de gènes hors fig1 ciblés par les FT ERK 
count_fig1_ERK_TF <- TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes %>%
  dplyr::filter(is_fig1 == "yes") %>% # Extraction des gènes cibles, de la figure 1 
  pull(targeted_by_ERK_TF) %>%
  str_count(",") # Compte le nombre de FT ERK ciblant chaque gène
count_fig1_ERK_TF <- count_fig1_ERK_TF + 1
count_fig1_ERK_TF[is.na(count_fig1_ERK_TF)] <- 0

count_not_fig1_ERK_TF <- TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes %>%
  dplyr::filter(is_fig1 != "yes") %>%
  pull(targeted_by_ERK_TF) %>%
  str_count(",") # Compte le nombre de FT ERK ciblant chaque gène
count_not_fig1_ERK_TF <- count_not_fig1_ERK_TF + 1
count_not_fig1_ERK_TF[is.na(count_not_fig1_ERK_TF)] <- 0

p1 <- ggplot() +
  aes(count_fig1_ERK_TF) +
  geom_histogram(binwidth=1, colour="black", fill="white") +
  xlim(-1, 30) +
  xlab("Number of genes from fig1 targeted by ERK TFs")
p2 <- ggplot() +
  aes(count_not_fig1_ERK_TF) +
  geom_histogram(binwidth=1, colour="black", fill="white") +
  xlim(-1, 30) + 
  xlab("Number of genes not from fig1 targeted by ERK TFs")

pdf("count_ERK_TF_targets.pdf")
p <- grid.arrange(p1, p2, nrow = 2)
p
dev.off()

## Heatmaps
# Heatmap rerpésentant la proportion des gènes de la figure 1 qui sont ciblés par les FT ERK pour chaque voie de signalisation
# En colonne : les 4 voies de signalisation ; En ligne : les FT ERK

# Extraction des FT ERK
TF_ERK_genes <- TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes %>%
  dplyr::filter(is_ERK == "yes", TF2DNA_db == "yes") %>%
  pull(symbol)
names(TF_ERK_genes) <- TF_ERK_genes 

# Extraction des listes de gènes des 4 voies de signalisation (SYMBOL + ENSEMBL)
xlsx_pathways_list <- list("angiogenesis_upreg" = "Angiogenesis (upreg)",
                           "Inflammation_upreg" = "Inflammation (upreg)",
                           "Angiogenesis_downreg" = "Angiogenesis (downreg)",
                           "Cell_adhesion_downreg" = "Cell adhesion (downreg)")

get_gene_list_xlsx <- function(sheet_name){
  read.xlsx("ERK_pathway_db/LRRK2_all datasets compiled_FINAL.xlsx", sheet_name, colNames = FALSE) %>%
  unique %>%
  pull(X1) %>%
  mapIds(x = org.Hs.eg.db, keytype = "ENSEMBL", column = "SYMBOL") %>%
  as.data.frame %>%
  rownames_to_column("ENSEMBL") %>%
  `colnames<-` (c("ENSEMBL", "SYMBOL")) %>%
  arrange(SYMBOL)
}

xlsx_pathways_genes <- map(xlsx_pathways_list, ~ get_gene_list_xlsx(.x))

# Extraction des listes de FT ERK pour tous les gènes des 4 voies de signalisation
get_pathway_ERK_TF <- function(pathway_df){
  TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes %>%
    dplyr::filter(symbol %in% pathway_df$SYMBOL) %>%
    pull(targeted_by_ERK_TF)
}
pathways_ERK_TF_list <- map(xlsx_pathways_genes, ~ get_pathway_ERK_TF(.x))

# Calcul de la proportion de gènes ciblés par les FT ERK dans chaque voie de signalisation
get_TF_proportion <- function(pathways_ERK_TF, symbol){
  symbol_is_TF <- str_detect(pathways_ERK_TF, symbol)
  symbol_is_TF[is.na(symbol_is_TF)] <- FALSE
  nb_TF <- as.vector(table(symbol_is_TF)["TRUE"])
  if(is.na(nb_TF)){
    prop_TF = 0
  }else{
  prop_TF <- nb_TF / length(pathways_ERK_TF) * 100
  }
  prop_TF
}

get_TF_proportion_df <- function(TF_symbol){
  cbind("symbol" = TF_symbol, map_dfr(pathways_ERK_TF_list, ~ get_TF_proportion(.x, TF_symbol)))
}

TF_proportion_df <- map_dfr(TF_ERK_genes, ~ get_TF_proportion_df(.x)) %>%
  column_to_rownames("symbol") %>%
  as.matrix

pdf("heatmap_proportion_TF_ERK_in_pathways.pdf")
hm <- Heatmap(TF_proportion_df,
              column_names_rot = 45,
              name = "% Fig1 genes\ntargeted",
              column_title = "Proportion of Figure 1 genes targeted\nby ERK FTs for each signaling pathway")
hm <- draw(hm)
dev.off()

## Graphe d'interaction des FT ERK sur les gènes de la figure 1
# Générer l'input pour la table pour Cytoscape, cad une ligne par interaction
ERK_TF_targets <- read.xlsx("livrables/livrables_20220602/TF_target_ERK_Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes.xlsx",
                            sheet = "ERK_TF_targets") %>%
  dplyr::filter(target_name %in% fig1_genes$SYMBOL) %>%
  write_csv("cytoscape/ERK_fig1_edges_table.csv")

# Ajouter un fichier de metadonnées, décrivant si les gènes (noeuds) sont issus de la fig1, ou de ERK, ainsi que leur voie de signalisation
nodes_table <- TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes$symbol %>%
  tibble(symbol = .) %>%
  mutate(angiogenesis_upreg = ifelse(symbol %in% xlsx_pathways_genes$angiogenesis_upreg$SYMBOL, "angiogenesis_upreg", NA)) %>%
  mutate(Inflammation_upreg = ifelse(symbol %in% xlsx_pathways_genes$Inflammation_upreg$SYMBOL, "inflammation_upreg", NA)) %>%
  mutate(Angiogenesis_downreg = ifelse(symbol %in% xlsx_pathways_genes$Angiogenesis_downreg$SYMBOL, "angiogenesis_downreg", NA)) %>%
  mutate(Cell_adhesion_downreg = ifelse(symbol %in% xlsx_pathways_genes$Cell_adhesion_downreg$SYMBOL, "cell_adhesion_downreg", NA)) %>%
  pivot_longer(-symbol) %>%
  na.omit %>%
  group_by(symbol) %>%
  summarize(pathways = paste(value, collapse = ",")) %>%
  full_join(TF_target_list$Cortical_astrocytes_vs_Control_iPSC_derived_astrocytes) %>%
  dplyr::select(symbol, is_ERK, is_fig1, pathways) %>%
  dplyr::filter(is_ERK == "yes" | is_fig1 == "yes") %>%
  mutate(is_ERK_fig1 = ifelse(is_ERK == "yes" & is_fig1 != "yes", "ERK",
                       ifelse(is_ERK != "yes" & is_fig1 == "yes", "fig1", 
                       ifelse(is_ERK == "yes" & is_fig1 == "yes", "ERK_fig1", "none")))) %>%
  write_csv("cytoscape/fig1_gene_pathways.csv")

## Heatmap de la proportion des cibles FT ERK #1 qui sont également cibles du FT ERK #2
# Extraction des FT ERK qui ciblent les gènes de la fig1 (pas les autres FT)
ERK_TF_list <- unique(ERK_TF_targets$tf_name)

# Création de matrices contenant le nombre de gènes ciblés communs pour chaque FT ERK
compare_ERK_TF_targets_matrix <- matrix(NA, nrow=length(ERK_TF_list), ncol=length(ERK_TF_list), byrow=TRUE)
rownames(compare_ERK_TF_targets_matrix) <- ERK_TF_list
colnames(compare_ERK_TF_targets_matrix) <- ERK_TF_list

compare_ERK_TF_targets_matrix_percent <- compare_ERK_TF_targets_matrix

for (i in 1:length(ERK_TF_list)){
  tmp <- i + 1
  compare_ERK_TF_targets_matrix_percent[i,i] <- 100
  for (j in tmp:length(ERK_TF_list)){
    gene1 <- ERK_TF_list[i]
    gene2 <- ERK_TF_list[j]
    if(j != 30 & gene1 != gene2){
      list_target1 <- ERK_TF_targets$target_name[ERK_TF_targets$tf_name == gene1]
      list_target2 <- ERK_TF_targets$target_name[ERK_TF_targets$tf_name == gene2]
      compare_targets <- as.vector(table(list_target1 %in% list_target2)["TRUE"])
      compare_targets_percent <- compare_targets * 100 / nrow(fig1_genes)
      if(is.na(compare_targets)){
        compare_targets <- 0
        compare_targets_percent <- 0
      }
      compare_ERK_TF_targets_matrix[i,j] <- compare_targets
      compare_ERK_TF_targets_matrix[j,i] <- compare_targets
      compare_ERK_TF_targets_matrix_percent[i,j] <- compare_targets_percent
      compare_ERK_TF_targets_matrix_percent[j,i] <- compare_targets_percent
    }
  }
}

library(circlize)
col_fun = colorRamp2(c(0, 100), c("white", "darkblue"))

hm <- Heatmap(compare_ERK_TF_targets_matrix, column_names_rot = 45, col = c("white", "darkblue"))
pdf("hm_compare_ERK_TF_targets.pdf")
hm <- draw(hm)
dev.off()

hm <- Heatmap(compare_ERK_TF_targets_matrix_percent, column_names_rot = 45, col =  c("white", "darkblue"))
pdf("hm_compare_ERK_TF_targets_percent.pdf")
hm <- draw(hm)
dev.off()
