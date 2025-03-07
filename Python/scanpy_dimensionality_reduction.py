
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/dimensionality_reduction.html

import scanpy as sc

sc.settings.verbosity = 0  # 0, 1, 2 (default), 3
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

adata = sc.read(
    filename="./data/s4d8_quality_control_normalization_feature_selection.h5ad"
)

# Set the main expression matrix to log1p normalized values (Shifted Logarithm Normalization)
adata.layers["original_counts"] = adata.X.copy()  # Store raw count before overwriting
adata.X = adata.layers["log1p_norm"]  # Use log1p normalized counts for analysis

# Mark highly variable genes as highly deviant ones
# This allows the use of `use_highly_variable=True` in PCA
adata.var["highly_variable"] = adata.var["highly_deviant"]

# Method 1: PCA
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)

# Plot PCA results (scatter plot of PC1 vs PC2)
sc.pl.pca_scatter(adata, color="total_counts")

# Method 2: t-SNE
# Use PCA-transformed data as input to speed up t-SNE
sc.tl.tsne(adata, use_rep="X_pca")

# Plot t-SNE results
sc.pl.tsne(adata, color="total_counts")

# Method 3: UMAP
# Compute the neighborhood graph (required for UMAP)
sc.pp.neighbors(adata)

# UMAP
sc.tl.umap(adata)

# Plot UMAP results
sc.pl.umap(adata, color="total_counts")

# Plot UMAP with multiple quality control (QC) metrics
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
)

adata.write("./data/s4d8_quality_control_normalization_feature_selection_dimensionality_reduction.h5ad")
