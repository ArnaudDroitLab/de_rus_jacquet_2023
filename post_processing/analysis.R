setwd("/mnt/tn02_bioinfo/JACA/IPS/post_processing/")
.libPaths(c("/home/cjbparlant/R/vls137", .libPaths()))
library(tidyverse)
library(rnaseq)
library(openxlsx)
library(factoextra)
library(ggrepel)
library(ComplexHeatmap)

# Import
filenames <- rnaseq::get_filenames("../rnaseq_pipeline/results/kallisto/")
filenames <- filenames[which(!names(filenames) %in% c("SRR7406472", "SRR7406473", "SRR7406474", "SRR7406475"))]

anno <- "/is3/projects/PUBLIC/rnaseq_anno/org/Hs/Hs.Ensembl104.csv"

metadata <- read_csv("../metadata/SraRunTable_all.csv")

txi <- import_kallisto(filenames, anno = anno)
saveRDS(txi, "r_objects/txi.rds")
# txi <- readRDS("r_objects/txi.rds")

raw_counts <- get_anno_df(txi, "raw_counts")
tpm <- get_anno_df(txi, "tpm")

write_csv(raw_counts, "livrables/raw_counts.csv")
write.xlsx(raw_counts, "livrables/raw_counts.xlsx")

write_csv(tpm, "livrables/tpm.csv")
write.xlsx(tpm, "livrables/tpm.xlsx")

# ------------------------------------ PCA ----------------------------------------------
pdf("livrables/pca.pdf")
res_pca <- produce_pca(txi)
dev.off()

res_pca_df <- res_pca$df %>%
  left_join(metadata, by = c("sample" = "Run"))
res_pca_df$disease_state[is.na(res_pca_df$disease_state)] <- res_pca_df$Cell_type[is.na(res_pca_df$disease_state)]

p <- ggplot2::ggplot(res_pca_df, ggplot2::aes_string(x = "Dim1", y = "Dim2", color = "Cell_type", shape = "disease_state")) +
     ggplot2::geom_point(size = 3) +
     ggrepel::geom_text_repel(ggplot2::aes(label = sample), color = "black", force = 10) +
     ggplot2::theme_bw() +
     ggplot2::xlab(paste0("Dim1 (", res_pca$pca$eig[1, 2] %>% round(2), "%)")) +
     ggplot2::ylab(paste0("Dim2 (", res_pca$pca$eig[2, 2] %>% round(2), "%)"))

pdf("livrables/pca_colored.pdf")
print(p)
dev.off()

# ------------------------------------ HIERARCHICAL CLUSTERING ----------------------------------------------
res.hcpc <- FactoMineR::HCPC(res_pca$pca, graph = FALSE)
pdf("livrables/hierarchical_clustering.pdf")
fviz_dend(res.hcpc, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE,
           rect_border = "jco", labels_track_height = 0.8)
dev.off()

# ------------------------------------ DE ----------------------------------------------
design <- read_csv("design.csv")

dds <- deseq2_analysis(txi, design, ~ group)
saveRDS(dds, "r_objects/dds.rds")

res <- DESeq2::results(dds, contrast = c("group", "astrocytes", "iPSC")) %>%
    as.data.frame %>%
    rownames_to_column("id") %>%
    left_join(txi$anno, by = "id") %>%
    dplyr::select(-id) %>%
    dplyr::select(ensembl_gene:transcript_type, everything()) %>%
    arrange(padj)
write_csv(res, "livrables/de_astrocytes_vs_iPSC.csv")
write.xlsx(res, "livrables/de_astrocytes_vs_iPSC.xlsx")

# ------------------------------------ VOLCANO PLOTS ----------------------------------------------
p <- produce_volcano(res, fc_threshold = 1.5)$p +
    geom_text_repel(data = head(res, 10), mapping = aes(x = log2FoldChange, y = -log10(padj), label = symbol), color = "black")
pdf("livrables/volcano_astrocytes_vs_iPSC_FC_1_5.pdf")
print(p)
dev.off()

p <- produce_volcano(res, fc_threshold = 3)$p +
    geom_text_repel(data = head(res, 10), mapping = aes(x = log2FoldChange, y = -log10(padj), label = symbol), color = "black")
pdf("livrables/volcano_astrocytes_vs_iPSC_FC_3.pdf")
print(p)
dev.off()

# up_regulated <- res[res$log2FoldChange <= log2(1/3) & res$padj <= 0.05,] %>%
#     na.omit()
# write_csv(up_regulated, "livrables/up_regulated.csv")

# down_regulated <- res[res$log2FoldChange >= log2(3) & res$padj <= 0.05,] %>%
#     na.omit()
# write_csv(down_regulated, "livrables/down_regulated.csv")
