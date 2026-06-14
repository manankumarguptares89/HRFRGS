#!/usr/bin/env Rscript
#
# Phase 3: WGCNA Construction & Domain-Aware Pruning
# Builds weighted gene co-expression network and applies topology-based pruning
# Author: Biomarker Discovery Pipeline
# Date: 2025
#

library(WGCNA)
library(tidyverse)
library(igraph)
library(stringr)

# Allow multi-threading (optional but recommended)
allowWGCNAThreads()

# Set random seed for reproducibility
set.seed(42)
setwd("D:/Projects/Dissertation")

# Output directory
output_dir <- 'output'
dir.create(output_dir, showWarnings = FALSE)

cat("========================================\n")
cat("Phase 3: WGCNA Construction & Domain-Aware Pruning\n")
cat("========================================\n\n")

# ============================================================================
# HARD-CODED PARAMETERS (NO YAML)
# ============================================================================

# ---- WGCNA PARAMETERS ----
wgcna_config <- list(
  power = 0,                 # Soft thresholding power (auto-detect if 0)
  min_module_size = 20,      # Minimum genes per module
  deepSplit = 2,             # Sensitivity of module detection
  scale_free_r2 = 0.85       # Target scale-free topology fit
)

# ---- PRUNING PARAMETERS ----
pruning_config <- list(
  decay_lambda = 0.5,              # Exponential decay for PPI distance
  distance_threshold = 3,          # Max PPI distance allowed
  high_confidence_threshold = 0.95 # Always keep very strong correlations
)

# ============================================================================
# LOAD INPUT DATA
# ============================================================================
cat("Loading input data...\n")

# DEG results from Phase 1
degs_file <- file.path(output_dir, 'phase1_degs.csv')
if (!file.exists(degs_file)) {
  stop(sprintf("DEG file not found: %s", degs_file))
}

degs <- read.csv(degs_file, row.names = 1)
deg_genes <- rownames(degs[degs$Include_In_Network, ])
cat(sprintf("DEGs selected: %d\n", length(deg_genes)))

# Expression matrix
expr_file <- 'data/expression_matrix.csv'
expr_matrix <- read.csv(expr_file, row.names = 1)
expr_matrix <- expr_matrix[deg_genes, ]
cat(sprintf("Expression matrix: %d genes × %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# Filtered PPI network from Phase 2
ppi_file <- file.path(output_dir, 'phase2_filtered_ppi.csv')
if (!file.exists(ppi_file)) {
  warning("Filtered PPI not found — pruning will be skipped")
  ppi <- NULL
} else {
  ppi <- read.csv(ppi_file)
  cat(sprintf("Filtered PPI edges: %d\n", nrow(ppi)))
}

# ============================================================================
# STEP 1: NORMALIZATION
# ============================================================================
cat("\n--- Step 1: Data Normalization ---\n")

expr_vst <- log2(expr_matrix + 1)
cat("Applied log2 transformation\n")

# Remove zero-variance genes
expr_vst <- expr_vst[apply(expr_vst, 1, var) > 0, ]
cat(sprintf("Genes with variance > 0: %d\n", nrow(expr_vst)))

# ============================================================================
# STEP 2: CORRELATION MATRIX
# ============================================================================
cat("\n--- Step 2: Correlation Matrix ---\n")

correlation_matrix <- abs(cor(t(expr_vst), method = "pearson"))

cat(sprintf("Correlation matrix: %d × %d\n",
            nrow(correlation_matrix), ncol(correlation_matrix)))
cat(sprintf("Mean correlation: %.3f\n",
            mean(correlation_matrix[upper.tri(correlation_matrix)])))

# ============================================================================
# STEP 3: SOFT-THRESHOLDING POWER
# ============================================================================
cat("\n--- Step 3: Scale-Free Topology ---\n")

power_value <- wgcna_config$power
cat(sprintf("Using soft-thresholding power: %d\n", power_value))

# ============================================================================
# STEP 3: AUTOMATIC SOFT-THRESHOLDING POWER SELECTION
# ============================================================================
cat("\n--- Step 3: Scale-Free Topology (Auto-selection) ---\n")

# 1. Define a range of powers to test (typically 1 to 20 or 30)
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 2. Call the network topology analysis function
# Note: input needs to be the transposed expression matrix (Samples as rows)
sft <- pickSoftThreshold(
  t(expr_vst), 
  powerVector = powers, 
  verbose = 5, 
  networkType = "unsigned"
)

# 3. Select the power
# We look for the lowest power that hits the R^2 threshold (e.g., 0.80)
power_value <- sft$powerEstimate

# Fallback: If no power hits the threshold, use a default or the highest R^2
if (is.na(power_value)) {
  cat("Warning: Scale-free topology threshold not reached. Selecting power with highest R^2.\n")
  power_value <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
}

cat(sprintf("Selected optimal soft-thresholding power: %d\n", power_value))

# Optional: Save the SFT plot to check the "Elbow"
pdf(file.path(output_dir, "phase3_sft_plot.pdf"), width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=wgcna_config$scale_free_r2, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# ============================================================================
# STEP 4: ADJACENCY MATRIX
# ============================================================================
cat("\n--- Step 4: Adjacency Matrix ---\n")

adjacency_matrix <- correlation_matrix ^ power_value

cat(sprintf("Mean adjacency: %.4f\n", mean(adjacency_matrix)))

# ============================================================================
# STEP 5: TOM MATRIX
# ============================================================================
cat("\n--- Step 5: Topological Overlap Matrix (TOM) ---\n")

TOM <- TOMsimilarity(adjacency_matrix)
diag(TOM) <- 1

cat(sprintf("Mean TOM: %.4f\n", mean(TOM[upper.tri(TOM)])))

# ============================================================================
# STEP 6: MODULE DETECTION
# ============================================================================
cat("\n--- Step 6: Module Detection ---\n")

dissTOM <- 1 - TOM
hc <- hclust(as.dist(dissTOM), method = "average")

module_labels <- cutreeDynamic(
  dendro = hc,
  distM = dissTOM,
  deepSplit = wgcna_config$deepSplit,
  minClusterSize = wgcna_config$min_module_size,
  verbose = 0
)

gene_names <- rownames(expr_vst)

# Ensure length match
stopifnot(length(gene_names) == length(module_labels))

module_assignment <- data.frame(
  Gene = gene_names,
  Module = paste0("Module_", module_labels),
  stringsAsFactors = FALSE
)

rownames(module_assignment) <- gene_names


module_sizes <- table(module_assignment$Module)

cat(sprintf("Modules detected: %d\n", length(module_sizes)))
cat(sprintf("Module size range: %d–%d\n",
            min(module_sizes), max(module_sizes)))
# ============================================================================
# FASTER VECTORIZED STEP 7 (Recommended for 6k+ Genes)
# ============================================================================
cat("\n--- Step 7: Vectorized Domain-Aware Pruning ---\n")

if (!is.null(ppi)) {
  # 1. Sync TOM names
  gene_names <- rownames(expr_vst)
  rownames(TOM) <- gene_names
  colnames(TOM) <- gene_names
  #TOM <- adjacency_matrix
  
  # 2. Get distances (Matrix form)
  cat("Computing PPI distance matrix...\n")
  ppi_graph <- graph_from_data_frame(ppi[,1:2], directed=F)
  ppi_graph <- add_vertices(ppi_graph, nv=length(setdiff(gene_names, V(ppi_graph)$name)), 
                            name=setdiff(gene_names, V(ppi_graph)$name))
  ppi_dist <- distances(ppi_graph, v=gene_names, to=gene_names)
  
  # 3. Create Masking Matrices (Vectorized Math)
  cat("Applying pruning filters...\n")
  
  # A. Identify pairs to delete (Dist > Threshold)
  prune_mask <- ppi_dist > pruning_config$distance_threshold
  
  # B. Calculate Decay Matrix (Penalty for Dist > 1)
  decay_matrix <- matrix(1, nrow=nrow(TOM), ncol=ncol(TOM))
  decay_mask <- ppi_dist > 1 & ppi_dist <= pruning_config$distance_threshold
  decay_matrix[decay_mask] <- exp(-pruning_config$decay_lambda * (ppi_dist[decay_mask] - 1))
  
  # C. High-Confidence Protection Mask
  protect_mask <- correlation_matrix > pruning_config$high_confidence_threshold
  
  # 4. Final Calculation
  pruned_adjacency <- TOM * decay_matrix       # Apply decay
  pruned_adjacency[prune_mask] <- 0            # Cut distant edges
  pruned_adjacency[protect_mask] <- adjacency_matrix[protect_mask] # Restore high-confidence
  
  cat("Pruning complete.\n")
}

# ============================================================================
# STEP 8: SAVE RESULTS
# ============================================================================
cat("\n--- Step 8: Saving Results ---\n")

#save(TOM, file = file.path(output_dir, "phase3_tom_matrix.RData"))
write.csv(pruned_adjacency,
          file.path(output_dir, "phase3_pruned_adjacency.csv"),
          quote = FALSE)
write.csv(module_assignment,
          file.path(output_dir, "phase3_modules.csv"),
          quote = FALSE)

writeLines(
  c(
    sprintf("Soft Threshold Power: %d", power_value),
    sprintf("Modules detected: %d", length(module_sizes)),
    sprintf("Mean module size: %.1f", mean(module_sizes))
  ),
  file.path(output_dir, "phase3_power_info.txt")
)

cat("\n✓ Phase 3 completed successfully!\n")
