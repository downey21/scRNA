
# -*- coding: utf-8 -*-

rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

pbmc <- base::readRDS(file = "./data/pbmc3k_qc_feature_selection_dimensionality_reduction.rds")

# Leiden Clustering
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
pbmc <- Seurat::AddMetaData(pbmc, metadata = pbmc$seurat_clusters, col.name = "leiden_res0_25")

pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
pbmc <- Seurat::AddMetaData(pbmc, metadata = pbmc$seurat_clusters, col.name = "leiden_res0_5")

pbmc <- Seurat::FindClusters(pbmc, resolution = 1.0)
pbmc <- Seurat::AddMetaData(pbmc, metadata = pbmc$seurat_clusters, col.name = "leiden_res1")

p1 <- Seurat::DimPlot(pbmc, reduction = "umap", group.by = "leiden_res0_25", label = TRUE) + ggtitle("Leiden Clustering (res = 0.25)")
p2 <- Seurat::DimPlot(pbmc, reduction = "umap", group.by = "leiden_res0_5", label = TRUE) + ggtitle("Leiden Clustering (res = 0.5)")
p3 <- Seurat::DimPlot(pbmc, reduction = "umap", group.by = "leiden_res1", label = TRUE) + ggtitle("Leiden Clustering (res = 1.0)")

p1 + p2 + p3
