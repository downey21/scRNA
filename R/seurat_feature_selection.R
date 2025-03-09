
# -*- coding: utf-8 -*-

rm(list = ls())

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scry)
library(scran)
library(SingleCellExperiment)

pbmc <- base::readRDS(file = "./data/pbmc3k_qc_normalization.rds")

# Feature Selection
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Deviance Feature Selection

# Higher deviance:
# The gene exhibits a structured (non-random) expression pattern with potential biological significance.
# It likely shows meaningful differences across cell populations.
# It is less likely to be noise and more likely to contain important biological information.

# Lower deviance:
# The gene expression is relatively stable and predictable.
# It shows low variability and may not differentiate between cell populations.
# It could be a candidate for removal in downstream analysis.

sce <- Seurat::as.SingleCellExperiment(pbmc)

sce <- scry::devianceFeatureSelection(sce, assay = "counts")
binomial_deviance <- SingleCellExperiment::rowData(sce)$binomial_deviance

top_features <- names(sort(binomial_deviance, decreasing = TRUE)[1:4000])

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = binomial_deviance,
    col.name = "binomial_deviance"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = setNames(rownames(pbmc) %in% top_features, rownames(pbmc)),
    col.name = "highly_deviant"
)

head(pbmc@assays$RNA@meta.data)

# Scran Highly Variable Gene Selection
gene_var <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(gene_var, n = 2000)

var_data <- data.frame(
    mean = gene_var$mean,
    variance = gene_var$total,
    dispersion = gene_var$bio,
    highly_variable = rownames(sce) %in% hvgs
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = setNames(var_data$highly_variable, rownames(pbmc)),
    col.name = "highly_variable"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = setNames(var_data$mean, rownames(pbmc)),
    col.name = "means"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = setNames(var_data$dispersion, rownames(pbmc)),
    col.name = "dispersions"
)

head(pbmc@assays$RNA@meta.data)

ggplot2::ggplot(var_data, aes(x = mean, y = dispersion, color = highly_variable)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Highly Variable Genes: Mean vs Dispersion")

base::saveRDS(pbmc, file = "./data/pbmc3k_qc_normalization_feature_selection.rds")
