# Differential expression analysis
# Project: DE analysis in mice gastric cells
#
# By Juan Jovel (use it at your own risk)
# (juan.jovel@ucalgary.ca)
# 
# Last revision: Oct. 13, 2025


library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(gplots)
library(pheatmap)
library(biomaRt)
library(EnhancedVolcano)
library(vegan)
library(UpSetR)

rm(list = ls())

setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/markUngrin/PolinaKvasova_project')
data     <- read.table("allSamples_counts.tsv", sep = '\t', header = T, row.names = 1, stringsAsFactors = T)
metadata <- read.table("metadata.tsv", sep = '\t', header = T, row.names = 1)

data_sorted <- data[, row.names(metadata)]

# Calculate row means for each gene across all samples
row_means <- rowMeans(data_sorted)

# Count how many samples per gene have counts >= 5
samples_above_5 <- rowSums(data_sorted >= 5)

# Filter: keep genes with average >= 5 in at least 3 samples
data_filtered <- data_sorted[samples_above_5 >= 3, ]

# Convert to matrix
data_matrix <- as.matrix(round(data_filtered), digits = 0) 

# Import data matrix into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data_matrix, colData = metadata, design =~ group)

# Perform variance stabilizing transformation for PCA
vsd <- vst(dds, blind = FALSE)

# Generate PCA data
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Custom PCA plot with larger dots, black stroke, and custom fill colors
ggplot(pcaData, aes(x = PC1, y = PC2, fill = group, label = name)) +
  geom_point(size = 8, shape = 21, color = "black", stroke = 1.5) +
  geom_text(hjust = 0.5, vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Gastric Organoid Culture Conditions") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 11)
  ) +
  scale_fill_manual(values = c("01_primary" = "orchid",
                               "02_2D_smallMol" = "#3498DB", 
                               "03_3D_smallMol" = "firebrick",
                               "04_3D_purProt" = "forestgreen"))

# Save the plot
ggsave("PCA_by_group.pdf", width = 10, height = 7)

# Run DE analysis
dds <-  DESeq(dds, fitType='local')

# Distance between groups
# Use more than just PC1 and PC2 for distance calculation
pca_obj <- prcomp(t(assay(vsd)))
pca_scores <- as.data.frame(pca_obj$x[, 1:5])  # Use first 5 PCs
pca_scores$group <- metadata$group

# Calculate centroids using multiple PCs
centroids_multi <- aggregate(. ~ group, data = pca_scores, FUN = mean)
rownames(centroids_multi) <- centroids_multi$group
centroids_multi$group <- NULL

# Calculate distances using multiple PCs
dist_multi <- dist(centroids_multi, method = "euclidean")
dist_multi_df <- as.matrix(dist_multi)

print("Distances using first 5 PCs:")
print(round(dist_multi_df, 2))

# Create distance heatmap
pheatmap(dist_multi_df,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "white",
         color = colorRampPalette(c("white", "orchid", "brown"))(50),
         fontsize_number = 14,
         cellwidth = 60,
         cellheight = 60)

# Quantify distances
# Perform PCA
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Calculate centroids
centroids <- aggregate(cbind(PC1, PC2) ~ group, data = pcaData, FUN = mean)
rownames(centroids) <- centroids$group

# Calculate distances
dist_matrix <- dist(centroids[, c("PC1", "PC2")], method = "euclidean")
dist_df <- as.matrix(dist_matrix)

# Print results
cat("\n=== PCA Distance Analysis ===\n")
cat("\nPairwise Euclidean distances between group centroids (PC1 + PC2):\n")
print(round(dist_df, 2))

cat("\nDistances from Primary cells:\n")
primary_dist <- dist_df["01_primary", ]
primary_dist_sorted <- sort(primary_dist[primary_dist > 0])
print(round(primary_dist_sorted, 2))

cat("\nRanking of culture conditions by similarity to primary:\n")
cat("(Lower distance = more similar)\n")
for(i in 1:length(primary_dist_sorted)) {
  cat(sprintf("%d. %s: %.2f\n", i, names(primary_dist_sorted)[i], primary_dist_sorted[i]))
}

# Add centroids and distance lines to PCA plot
ggplot(pcaData, aes(x = PC1, y = PC2, fill = group)) +
  geom_point(size = 8, shape = 21, color = "black", stroke = 1.5, alpha = 0.7) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, fill = group), 
             size = 6, shape = 23, color = "blue", stroke = 2) +
  geom_segment(data = centroids[centroids$group != "01_primary", ],
               aes(x = centroids["01_primary", "PC1"], 
                   y = centroids["01_primary", "PC2"],
                   xend = PC1, yend = PC2),
               linetype = "dashed", color = "gray40", linewidth = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  scale_fill_manual(values = c("01_primary" = "orchid",
                               "02_2D_smallMol" = "#3498DB", 
                               "03_3D_smallMol" = "firebrick",
                               "04_3D_purProt" = "forestgreen"))

ggsave("PCA_with_centroids_and_distances.pdf", width = 12, height = 8)

# Extract the three comparisons of interest
# Use lfcShrink with the correct coefficient names
res_primary_vs_2D <- lfcShrink(dds, 
                               coef = "group_02_2D_smallMol_vs_01_primary", 
                               type = "apeglm")

res_primary_vs_3Dsm <- lfcShrink(dds, 
                                 coef = "group_03_3D_smallMol_vs_01_primary", 
                                 type = "apeglm")

res_primary_vs_3Dpp <- lfcShrink(dds, 
                                 coef = "group_04_3D_purProt_vs_01_primary", 
                                 type = "apeglm")

# View summaries
summary(res_primary_vs_2D)
summary(res_primary_vs_3Dsm)
summary(res_primary_vs_3Dpp)

# Get significant genes (adjusted p-value < 0.05)
sig_primary_vs_2D <- subset(res_primary_vs_2D, padj < 0.05 & abs(log2FoldChange) >= 1)
sig_primary_vs_3Dsm <- subset(res_primary_vs_3Dsm, padj < 0.05 & abs(log2FoldChange) >= 1)
sig_primary_vs_3Dpp <- subset(res_primary_vs_3Dpp, padj < 0.05 & abs(log2FoldChange) >= 1)

# Count significant genes
cat("Significant DEGs:\n")
cat("Primary vs 2D:", nrow(sig_primary_vs_2D), "\n")
cat("Primary vs 3D small mol:", nrow(sig_primary_vs_3Dsm), "\n")
cat("Primary vs 3D pur prot:", nrow(sig_primary_vs_3Dpp), "\n")

# Optional: Save results to files
write.table(as.data.frame(sig_primary_vs_2D), 
            file = "DESeq2_primary_vs_2D_smallMol.tsv", sep = "\t", quote = F)
write.table(as.data.frame(sig_primary_vs_3Dsm), 
            file = "DESeq2_primary_vs_3D_smallMol.tsv", sep = "\t", quote = F)
write.table(as.data.frame(sig_primary_vs_3Dpp), 
            file = "DESeq2_primary_vs_3D_purProt.tsv", sep = "\t", quote = F)

# Correlation analysis
expr <- assay(rld)
group_means <- aggregate(t(expr), by=list(metadata$group), FUN=mean)
rownames(group_means) <- group_means$Group.1
group_means <- group_means[,-1]
cor_matrix <- cor(t(group_means), method="pearson")

# Create heatmap with correlation coefficients displayed
pheatmap(cor_matrix, 
         main = "",
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "white",
         fontsize_number = 12,
         cellwidth = 60,
         cellheight = 60,
         color = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100))

###############################################
# UPSET PLOT - REPLACES VENN DIAGRAM
###############################################

# Get gene lists from significant results
genes_2D <- rownames(sig_primary_vs_2D)
genes_3Dsm <- rownames(sig_primary_vs_3Dsm)
genes_3Dpp <- rownames(sig_primary_vs_3Dpp)

# Create a list of all unique genes
all_genes <- unique(c(genes_2D, genes_3Dsm, genes_3Dpp))

# Create binary matrix for UpSet plot
upset_data <- data.frame(
  gene = all_genes,
  Primary_vs_2D = as.integer(all_genes %in% genes_2D),
  Primary_vs_3DsmallMol = as.integer(all_genes %in% genes_3Dsm),
  Primary_vs_3DpurProt = as.integer(all_genes %in% genes_3Dpp)
)

# Generate UpSet plot
png("DE_upset_plot.png", width = 800, height = 600)
upset(upset_data, 
      sets = c("Primary_vs_2D", "Primary_vs_3DsmallMol", "Primary_vs_3DpurProt"),
      order.by = "freq",
      keep.order = TRUE,
      main.bar.color = "steelblue",
      sets.bar.color = "darkgreen",
      text.scale = c(1.5, 1.3, 1.3, 1.3, 1.5, 1.3),
      point.size = 3.5,
      line.size = 1)
dev.off()

# Print overlap statistics
cat("\n=== Gene Overlap Statistics ===\n")
cat("Genes in 2D only:", sum(upset_data$Primary_vs_2D == 1 & 
                               upset_data$Primary_vs_3DsmallMol == 0 & 
                               upset_data$Primary_vs_3DpurProt == 0), "\n")
cat("Genes in 3D small mol only:", sum(upset_data$Primary_vs_2D == 0 & 
                                         upset_data$Primary_vs_3DsmallMol == 1 & 
                                         upset_data$Primary_vs_3DpurProt == 0), "\n")
cat("Genes in 3D pur prot only:", sum(upset_data$Primary_vs_2D == 0 & 
                                        upset_data$Primary_vs_3DsmallMol == 0 & 
                                        upset_data$Primary_vs_3DpurProt == 1), "\n")
cat("Genes shared by all three:", sum(upset_data$Primary_vs_2D == 1 & 
                                        upset_data$Primary_vs_3DsmallMol == 1 & 
                                        upset_data$Primary_vs_3DpurProt == 1), "\n")

###############################################

dist_mat <- dist(t(expr))
hc <- hclust(dist_mat)
plot(hc, labels=metadata$label, main="Hierarchical Clustering of Transcriptomes")

primary_expr <- rowMeans(expr[, metadata$group == "01_primary"])
cor(rank(primary_expr), rank(group_means["02_2D_smallMol", ]), method="spearman")

primary_expr <- rowMeans(expr[, metadata$group == "01_primary"])
cor(rank(primary_expr), rank(group_means["03_3D_smallMol", ]), method="spearman")

primary_expr <- rowMeans(expr[, metadata$group == "01_primary"])
cor(rank(primary_expr), rank(group_means["04_3D_purProt", ]), method="spearman")

##############################################
#           STARTS MACHINE LEARNING          #
##############################################
dds <- DESeqDataSetFromMatrix(countData = data_matrix, colData = metadata, design =~ group)
vst_counts <- assay(vst(dds, blind=TRUE))
vst_counts <- vst_counts[apply(vst_counts, 1, var) > 0.5, ]

labels <- metadata$group
X <- t(vst_counts)
y <- factor(make.names(metadata$group))

library(caret)
set.seed(123)

train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 10, 
                              classProbs = TRUE, summaryFunction = multiClassSummary)

rf_model <- train(
  x = X,
  y = y,
  method = "rf",
  trControl = train_control,
  tuneLength = 5,
  importance = TRUE
)

rf_model

###
pred <- predict(rf_model, X, type = "prob")
aggregate(pred$X01_primary, by=list(group=y), mean)

# Extract feature importance
importance <- varImp(rf_model)
plot(importance, top=30)

library(umap)
umap_res <- umap(X, n_neighbors = 5, random_state = 123)

plot(umap_res$layout, col=as.factor(y), pch=19)
###

pred <- predict(rf_model, X, type = "prob")
pred_df <- data.frame(sample = rownames(X), metadata, pred)

# Mean probability per group
aggregate(pred_df$X01_primary, by=list(pred_df$group), mean)

pred_class <- predict(rf_model, X)
confusionMatrix(pred_class, y)

library(umap)
umap_res <- umap(X, n_neighbors = 5, random_state = 123)
plot(umap_res$layout, col=as.factor(y), pch=19, cex=2,
     main="UMAP Projection of Gastric Organoid Conditions")
text(umap_res$layout, labels=metadata$label, pos=3)

importance <- varImp(rf_model)
plot(importance, top=30)

library(dplyr)
similarity_scores <- pred_df %>%
  group_by(group) %>%
  summarise(
    mean_primary_prob = mean(X01_primary),
    sd_primary_prob = sd(X01_primary)
  ) %>%
  arrange(desc(mean_primary_prob))
similarity_scores

library(stats)
coords <- data.frame(umap_res$layout, group = y)
centroids <- aggregate(coords[,1:2], by=list(coords$group), mean)
rownames(centroids) <- centroids$Group.1

dist_to_primary <- sqrt(
  (centroids["X01_primary",2] - centroids[,2])^2 +
    (centroids["X01_primary",3] - centroids[,3])^2
)

plot(varImp(rf_model), top = 20, layout = c(1, 4))

imp <- varImp(rf_model)$importance
imp$gene <- rownames(imp)

imp_long <- imp %>%
  pivot_longer(starts_with("X"), names_to = "group", values_to = "importance")

# top N
topN <- 20
imp_top <- imp_long %>%
  group_by(group) %>%
  slice_max(order_by = importance, n = topN, with_ties = FALSE) %>%
  ungroup()

# facet labels + order
lab <- c(
  X01_primary = "Primary",
  X02_2D_smallMol = "2D smallMol",
  X03_3D_smallMol = "3D smallMol",
  X04_3D_purProt  = "3D purProt"
)
imp_top$group <- factor(imp_top$group, levels = names(lab), labels = lab)

# order genes within facet
imp_top <- imp_top %>%
  group_by(group) %>% mutate(gene = reorder(gene, importance)) %>% ungroup()

# lollipop base
p <- ggplot(imp_top, aes(x = gene, y = importance)) +
  geom_segment(aes(xend = gene, y = 0, yend = importance),
               linewidth = 0.5, color = "grey50") +
  geom_point(size = 2.6, color = "#1e90ff") +
  facet_wrap(~ group, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Importance",
       title = "Random Forest Gene Importance per Group (top 20/class)") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 10, 10, 80)
  ) +
  scale_y_continuous(limits = c(-5, 100), expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = "off")

# add labels ONLY for Primary, placed in the small negative zone
p + geom_text(
  data = dplyr::filter(imp_top, group == "Primary"),
  aes(label = gene, y = -1.5),
  hjust = 1, size = 3
)

##############################################
#           ENDS MACHINE LEARNING            #
##############################################

##### ANNOTATION OF DE FEATURES #####
library(biomaRt)

# Create remote connection
my_mart <- useEnsembl('ensembl', dataset = "mmusculus_gene_ensembl")

# Make list of attributes to retrieve
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "entrezgene_id",
                "external_gene_name",
                "wikigene_description",
                "name_1006",
                "definition_1006",
                "namespace_1003")

# Function to extract records
pullRecords <- function(attributes, mart, filter_values){
  records <- getBM(attributes = attributes, filters = "ensembl_transcript_id", 
                   values = filter_values, mart = mart)
  return(records)
}

# Extract first hit only
getFirstMatch <- function(records, transcripts) {
  first_annotation_df <- data.frame()
  for (transcript in transcripts) {
    hits <- which(records$ensembl_transcript_id == transcript)
    if (length(hits) > 0) {
      first <- hits[1]
      first_annotation <- records[first,]
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    } else {
      first_annotation <- c(transcript, "_", "-", "-", "-", "-", "-", "-")
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    }
  }
  colnames(first_annotation_df) <- c("ensembl_transcript_id", "ensembl_gene_id", 
                                     "entrezgene_id", "external_gene_name", 
                                     "wikigene_description", "GO_group", 
                                     "GO_definition", "Ontology")
  return(first_annotation_df)
}

makeVolcanoPlot <- function(df, vp_file, title_text){
  keyvals <- ifelse(
    df$log2FoldChange < -1, 'forestgreen',
    ifelse(df$log2FoldChange > 1, 'firebrick1',
           'dodgerblue1'))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'firebrick1'] <- 'High'
  names(keyvals)[keyvals == 'dodgerblue1'] <- 'small FC'
  names(keyvals)[keyvals == 'forestgreen'] <- 'Low'
  
  evp <- EnhancedVolcano(df,
                         lab = df$external_gene_name,
                         x = 'log2FoldChange',
                         y = 'padj',
                         pCutoff = 0.05,
                         FCcutoff = 1,
                         colCustom = keyvals,
                         title = title_text,
                         subtitle = NULL,
                         colAlpha = 0.85,
                         shape = 20,
                         pointSize = 2,
                         labSize = 3)
  
  png(vp_file, width = 800, height = 600)
  print(evp)
  dev.off()
}

###############################################
# PROCESS ALL THREE COMPARISONS
###############################################

# Extract transcript IDs and clean them
transcripts1 <- row.names(sig_primary_vs_2D)
transcripts2 <- row.names(sig_primary_vs_3Dsm)
transcripts3 <- row.names(sig_primary_vs_3Dpp)

transcripts_clean1 <- gsub("\\.\\d+$", "", transcripts1)
transcripts_clean2 <- gsub("\\.\\d+$", "", transcripts2)
transcripts_clean3 <- gsub("\\.\\d+$", "", transcripts3)

# Pull annotations for each comparison
cat("\n=== Pulling annotations from Ensembl ===\n")
cat("Processing Primary vs 2D...\n")
transc_annotations1 <- pullRecords(attributes, my_mart, transcripts_clean1)
annotation_df1 <- getFirstMatch(transc_annotations1, transcripts_clean1)

cat("Processing Primary vs 3D small mol...\n")
transc_annotations2 <- pullRecords(attributes, my_mart, transcripts_clean2)
annotation_df2 <- getFirstMatch(transc_annotations2, transcripts_clean2)

cat("Processing Primary vs 3D pur prot...\n")
transc_annotations3 <- pullRecords(attributes, my_mart, transcripts_clean3)
annotation_df3 <- getFirstMatch(transc_annotations3, transcripts_clean3)

# Create list of results and annotations for iteration
results_list <- list(
  res_primary_vs_2D = list(result = res_primary_vs_2D, 
                           sig_result = sig_primary_vs_2D,
                           annotation = annotation_df1),
  res_primary_vs_3Dsm = list(result = res_primary_vs_3Dsm, 
                             sig_result = sig_primary_vs_3Dsm,
                             annotation = annotation_df2),
  res_primary_vs_3Dpp = list(result = res_primary_vs_3Dpp, 
                             sig_result = sig_primary_vs_3Dpp,
                             annotation = annotation_df3)
)

# Process each comparison
cat("\n=== Processing and saving annotated results ===\n")
for (res_name in names(results_list)) {
  prefix <- gsub("res_", "", res_name)
  cat(paste("Processing", prefix, "...\n"))
  
  result <- results_list[[res_name]]$result
  sig_result <- results_list[[res_name]]$sig_result
  annotation_df <- results_list[[res_name]]$annotation
  
  # Combine full results with annotations for volcano plot
  annotated_results <- cbind(as.data.frame(result), annotation_df)
  
  # Combine significant results with annotations
  annotated_sig_results <- cbind(as.data.frame(sig_result), annotation_df)
  annotated_sig_results <- annotated_sig_results[order(-annotated_sig_results$log2FoldChange), ]
  annotated_sig_results <- cbind(transcript=rownames(annotated_sig_results), 
                                 annotated_sig_results)
  
  # Save annotated significant results
  file_name <- paste(prefix, "q0.05_FC1_annotated.tsv", sep = '_')
  write.table(annotated_sig_results, file_name, quote = F, sep = '\t', row.names = F)
  
  # Create volcano plot
  vp_file <- paste(prefix, "volcanoPlot.png", sep = '_')
  makeVolcanoPlot(annotated_results, vp_file, 
                  title_text = paste("Volcano Plot:", gsub("_", " ", prefix)))
  
  cat(paste("Number of significant transcripts deregulated:", 
            nrow(annotated_sig_results), "\n"))
  cat(paste("Saved:", file_name, "\n"))
  cat(paste("Saved:", vp_file, "\n\n"))
}

cat("\n=== Analysis Complete ===\n")