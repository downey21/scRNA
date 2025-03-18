
# -*- coding: utf-8 -*-

# Normalization
# Normalization adjusts for differences in count depth caused by sequencing variability
# Normalization step ensures that statistical methods can be more accurately applied in subsequent analyses.

rm(list = ls())

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scran)
library(SingleCellExperiment)
library(BiocParallel)

pbmc <- readRDS(file = "./data/pbmc3k_qc.rds")

# Plot raw total counts distribution
# ggplot2::ggplot(pbmc@meta.data, aes(x = nCount_RNA)) +
#     ggplot2::geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
#     ggplot2::theme_minimal() +
#     ggplot2::ggtitle("Total Counts Distribution")

# Shifted Logarithm Normalization
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

raw_counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
log1p_norm_counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "data")

normalized_nCount_RNA <- Matrix::colSums(log1p_norm_counts)
pbmc[["log1p_norm"]] <- normalized_nCount_RNA[rownames(pbmc@meta.data)]

# Scran’s Pooling-Based Size Factor Estimation Method

# Scran is a normalization method that estimates size factors to account for differences in cell size.
# It better adjusts for cell-to-cell variation by computing size factors specific to each cell,
# allowing for more effective correction of technical variation.

# Explanation:
# - Shifted Logarithm: Applies a log1p transformation to stabilize variance in counts.
# - Scran Normalization: Uses a clustering-based method to estimate size factors, 
#   which adjusts for differences in cell-specific biases before applying log1p normalization.

# In most cases, Shifted Logarithm Normalization is sufficient.

sce <- Seurat::as.SingleCellExperiment(pbmc)

pbmc_pp <- pbmc
pbmc_pp <- Seurat::FindVariableFeatures(pbmc_pp, selection.method = "vst", nfeatures = 2000)
pbmc_pp <- Seurat::ScaleData(pbmc_pp)
pbmc_pp <- Seurat::RunPCA(pbmc_pp, npcs = 15)
pbmc_pp <- Seurat::FindNeighbors(pbmc_pp, dims = 1:15)
pbmc_pp <- Seurat::FindClusters(pbmc_pp, resolution = 0.5)
sce$clusters <- pbmc_pp$seurat_clusters

sce <- scran::computeSumFactors(
    sce,
    clusters = sce$clusters,
    min.mean = 0.1,
    BPPARAM = BiocParallel::MulticoreParam()
)

size_factors <- SingleCellExperiment::sizeFactors(sce)

pbmc@meta.data$size_factors <- size_factors

scran_norm_counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts") / size_factors
scran_norm_counts <- log1p(scran_norm_counts)

pbmc <- Seurat::SetAssayData(pbmc, assay = "RNA", layer = "data_scran", new.data = scran_norm_counts)

scran_normalized_nCount_RNA <- Matrix::colSums(scran_norm_counts)
pbmc[["scran_normalization"]] <- scran_normalized_nCount_RNA[rownames(pbmc@meta.data)]

p1 <- ggplot2::ggplot(pbmc@meta.data, aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Total counts (Raw)")

p2 <- ggplot2::ggplot(pbmc@meta.data, aes(x = log1p_norm)) +
    ggplot2::geom_histogram(bins = 100, fill = "tomato", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Shifted Logarithm Normalization")

p3 <- ggplot2::ggplot(pbmc@meta.data, aes(x = scran_normalization)) +
    ggplot2::geom_histogram(bins = 100, fill = "green", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Scran Normalization")

p1 + p2 + p3

base::saveRDS(pbmc, file = "./data/pbmc3k_qc_normalization.rds")
