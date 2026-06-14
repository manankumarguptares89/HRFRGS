#!/usr/bin/env Rscript
#
# Phase 1: Differential Expression Analysis using DESeq2
# Identifies Differentially Expressed Genes (DEGs) between Tumor and Normal samples
# Author: Biomarker Discovery Pipeline
# Date: 2025
#

library(DESeq2)
library(tidyverse)
library(stringr)
library(yaml)
library(dplyr)

setwd("D:/Projects/Dissertation")

# Configuration
config_file <- 'config/pipeline_config.yaml'
output_dir <- 'output'

# Create output directory
#dir.create(output_dir, showWarnings = FALSE)

# Load configuration
load_config <- function(file) {
  config <- yaml::read_yaml(file)
  return(config)
}

config <- load_config(config_file)
deseq_config <- config$deseq2

# Load input data
cat("Loading input data...\n")

# Expression matrix (genes x samples)
expr_file <- 'data/expression_matrix.csv'
if (!file.exists(expr_file)) {
  stop(sprintf("Expression file not found: %s", expr_file))
}

expr_matrix <- read.csv(expr_file, row.names = 1)
cat(sprintf("Expression matrix: %d genes x %d samples\n", 
            nrow(expr_matrix), ncol(expr_matrix)))

# Metadata (sample annotations)
meta_file <- 'data/metadata.csv'
if (!file.exists(meta_file)) {
  stop(sprintf("Metadata file not found: %s", meta_file))
}

metadata <- read.csv(meta_file, row.names = 1,check.names = FALSE)
cat(sprintf("Metadata: %d samples\n", nrow(metadata)))
# Convert expression sample names: dots → dashes
colnames(expr_matrix) <- gsub("\\.", "-", colnames(expr_matrix))

# Verify sample names match
if (!all(colnames(expr_matrix) %in% rownames(metadata))) {
  stop("Sample names in expression matrix don't match metadata")
}

# Reorder metadata to match expression matrix
metadata <- metadata[colnames(expr_matrix), ,drop= FALSE ]

#colnames(metadata) <- "Condition"
#metadata$Condition <- ifelse(metadata$Condition == 1, "Tumor", "Normal")
metadata$Condition <- factor(metadata$Status, levels = c("Normal", "Tumor")
)

stopifnot(
  is.data.frame(metadata),
  nrow(metadata) == ncol(expr_matrix),
  all(colnames(expr_matrix) == rownames(metadata)),
  sum(is.na(metadata$Condition)) == 0
)

cat("Metadata successfully converted: Normal/Tumor.\n")

cat("Metadata is valid for DESeq2.\n")
# Create DESeq2 object
cat("Creating DESeq2 object...\n")

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(expr_matrix),
  colData = metadata,
  design = ~ Condition
)

# Set reference level (Normal should be reference)
dds$Condition <- relevel(dds$Condition, ref = "Normal")

# Pre-filtering: remove genes with very low counts
cat("Pre-filtering low-abundance genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Genes retained after filtering: %d\n", nrow(dds)))

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get results
results <- results(dds, contrast = c("Condition", "Tumor", "Normal"))

# Convert to dataframe
results_df <- as.data.frame(results)
results_df$Gene <- rownames(results_df)
results_df <- results_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", 
                              "stat", "pvalue", "padj")]

# Apply filtering criteria
padj_threshold <- deseq_config$padj_threshold
lfc_threshold <- deseq_config$log2fc_threshold
base_mean_threshold <- deseq_config$base_mean_threshold

cat(sprintf("Filtering criteria:\n"))
cat(sprintf("  padj < %f\n", padj_threshold))
cat(sprintf("  |log2FC| > %f\n", lfc_threshold))
cat(sprintf("  base mean > %f\n", base_mean_threshold))

results_df <- results_df %>%
  mutate(
    Include_In_Network = (!is.na(padj)) & 
                         (padj < padj_threshold) & 
                         (abs(log2FoldChange) > lfc_threshold) &
                         (baseMean > base_mean_threshold)
  ) %>%
  arrange(padj)

n_degs <- sum(results_df$Include_In_Network)
n_up <- sum(results_df$Include_In_Network & results_df$log2FoldChange > 0)
n_down <- sum(results_df$Include_In_Network & results_df$log2FoldChange < 0)

cat(sprintf("\nDifferential Expression Results:\n"))
cat(sprintf("  Total DEGs: %d\n", n_degs))
cat(sprintf("  Upregulated: %d\n", n_up))
cat(sprintf("  Downregulated: %d\n", n_down))

# Show top DEGs
cat("\nTop 10 DEGs (by adjusted p-value):\n")
top_degs <- results_df %>% 
  filter(Include_In_Network) %>% 
  head(10) %>%
  dplyr::select(Gene, log2FoldChange, padj)

print(knitr::kable(top_degs), quote = FALSE)

# Save results
output_file <- file.path(output_dir, 'phase1_degs.csv')
write.csv(results_df, output_file, row.names = FALSE, quote = FALSE)
cat(sprintf("\nResults saved: %s\n", output_file))

# Generate volcano plot
cat("Generating volcano plot...\n")

results_df <- results_df %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-300),
    color = case_when(
      !Include_In_Network ~ "gray",
      log2FoldChange > 0 ~ "red",
      TRUE ~ "blue"
    )
  )

png(file.path(output_dir, 'phase1_volcano_plot.png'), 
    width = 800, height = 600, res = 100)

plot(results_df$log2FoldChange, results_df$neg_log10_padj,
     main = "Volcano Plot: Tumor vs Normal",
     xlab = "log2(Fold Change)",
     ylab = "-log10(padj)",
     col = results_df$color,
     pch = 16,
     cex = 0.5)

abline(h = -log10(padj_threshold), col = "black", lty = 2, lwd = 1)
abline(v = c(-lfc_threshold, lfc_threshold), col = "black", lty = 2, lwd = 1)

legend("topright", 
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "gray"),
       pch = 16)

dev.off()
cat(sprintf("Volcano plot saved: %s\n", 
            file.path(output_dir, 'phase1_volcano_plot.png')))

cat("\nPhase 1 (DESeq2) completed successfully!\n")
