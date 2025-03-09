
# -*- coding: utf-8 -*-

rm(list = ls())

library(Seurat)
library(dplyr)

pbmc <- base::readRDS(file = "./data/pbmc3k_qc_normalization_feature_selection.rds")

# Scaling
pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))

# Dimensionality Reduction

# Method 1: PCA
pbmc <- Seurat::RunPCA(pbmc, features = VariableFeatures(pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
Seurat::VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
Seurat::DimPlot(pbmc, reduction = "pca") + NoLegend()
Seurat::DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# Seurat::DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# Seurat::ElbowPlot(pbmc)

# Method 2: t-SNE
pbmc <- Seurat::RunTSNE(pbmc, dims = 1:10)

# Method 3: UMAP
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)

# Plot with multiple quality control (QC) metrics
p1 <- Seurat::DimPlot(pbmc, reduction = "tsne", label = TRUE)
p2 <- Seurat::DimPlot(pbmc, reduction = "umap", label = TRUE)
p3 <- Seurat::FeaturePlot(pbmc, features = c("nCount_RNA", "percent.mt"), reduction = "umap")

p1 + p2 + p3

base::saveRDS(pbmc, file = "./data/pbmc3k_qc_normalization_feature_selection_dimensionality_reduction.rds")
