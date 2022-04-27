setwd("/mnt/tn02_bioinfo/JACA/IPS/post_processing/")
.libPaths(c("/home/cjbparlant/R/vls137", .libPaths()))
library(tidyverse)
library(rnaseq)
library(openxlsx)
library(factoextra)
library(ggrepel)
library(sva)

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

regulated_files <- dir(path = "livrables", pattern = "regulated", full.names = TRUE)
regulated_files <- regulated_files[!str_detect(regulated_files, "TF")]
names(regulated_files) <- str_remove(regulated_files, "livrables/") %>%
  str_remove(".xlsx")

regulated_list <- map(regulated_files, ~ read.xlsx(.x))

# DB 1 : The human Trancripition Factors db
TheHumanTrancriptionFactors_db <- read.csv("TheHumanTrancriptionFactors_db.csv") %>%
  dplyr::select(Ensembl.ID)

# DB 2 : Human TFDB
HumanTFDB_db <- read_tsv("HumanTFDB_db.txt") %>%
  dplyr::select(Ensembl)

# DB 3 : TF2DNA
pscan_files <- dir(path = "pscan_files", pattern = "pscan", full.names = TRUE, recursive = TRUE)

concat_all_pscan <- function(pscan_file, i, len){
  print(paste(i, "/", len))
  pscan_db <- str_remove(pscan_file, "pscan_files/Homo-sapiens_") %>%
    str_remove("/.*$")

  pscan_df <- read_tsv(pscan_file, show_col_types = FALSE) %>%
    mutate(database = pscan_db) %>%
    dplyr::select(tf_name, target_name, binding_score, p_value, database) %>%
    transform(tf_name = as.character(tf_name))
}

# TF2DNA_db <- map_dfr(pscan_files, ~ concat_all_pscan(.x, which(pscan_files == .x), length(pscan_files)))
# TF2DNA_db <- TF2DNA_db %>% write_csv("TF2DNA_db.csv")
TF2DNA_db <- read_csv("TF2DNA_db.csv")

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

imap(de_res_TF_list, ~ get_TF_targets(.x, .y, regulated_list, TF2DNA_filtered))