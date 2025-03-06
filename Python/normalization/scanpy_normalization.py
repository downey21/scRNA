
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/normalization.html

# Normalization
# It is the process of adjusting the variability caused by sampling effects
# to bring the data variance within a consistent range.
# This ensures that statistical methods can be more accurately applied in subsequent analyses.

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.sparse import issparse, csr_matrix

import anndata2ri  # AnnData (.h5ad) -> SingleCellExperiment (.RData)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

sc.settings.verbosity = 0  # 0, 1, 2 (default), 3
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

pandas2ri.activate()
anndata2ri.activate()

adata = sc.read(
    filename="./data/s4d8_quality_control.h5ad",
    backup_url="https://figshare.com/ndownloader/files/40014331",
)

print(np.sum(adata.X, axis=1).A1[:10])
print(adata.obs["total_counts"][:10].values)

# Plot raw total counts distribution
fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=ax)
ax.set_title("Total counts")
plt.show()

# Shifted Logarithm Normalization

# Normalize total counts and log1p transform
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

scales_counts_1e4 = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
adata.layers["log1p_norm_1e4"] = sc.pp.log1p(scales_counts_1e4["X"], copy=True)

# Plot comparison before and after normalization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
sns.histplot(adata.layers["log1p_norm_1e4"].sum(1), bins=100, kde=False, ax=axes[2])
axes[2].set_title("Shifted logarithm (defualt)")
plt.show()

# Scran’s Pooling-Based Size Factor Estimation Method

# Scran is a normalization method that estimates size factors to account for differences in cell size.
# It better adjusts for cell-to-cell variation by computing size factors specific to each cell,
# allowing for more effective correction of technical variation.

# Preprocess data for Scran normalization
adata_pp = adata.copy()
print(np.sum(adata_pp.X, axis=1))

sc.pp.normalize_total(adata_pp)
print(np.sum(adata_pp.X, axis=1))

sc.pp.log1p(adata_pp)
print(np.sum(adata_pp.X, axis=1))

sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")

# Convert data matrix to appropriate format for R
data_mat = adata_pp.X.T
if issparse(data_mat):
    data_mat = data_mat.tocsc() if data_mat.nnz <= 2**31 - 1 else data_mat.tocoo()

# Pass data to R
data_mat_r = ro.r.matrix(data_mat.A if issparse(data_mat) else data_mat, 
                         nrow=data_mat.shape[0], 
                         ncol=data_mat.shape[1])
ro.globalenv["data_mat"] = data_mat_r
ro.globalenv["input_groups"] = ro.vectors.FactorVector(adata_pp.obs["groups"].tolist())

# Run Scran normalization in R
ro.r(
    """
    library(scran)
    library(SingleCellExperiment)
    library(BiocParallel)
    sce <- SingleCellExperiment(list(counts=data_mat))
    size_factors <- sizeFactors(
        computeSumFactors(
            sce,
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
        )
    )
    """
)

# Retrieve size factors from R to Python
size_factors = np.array(ro.r["size_factors"])
adata.obs["size_factors"] = size_factors

# Apply size factor normalization
scran_norm = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran_norm))

# Plot comparison before and after Scran normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()

fig, axes = plt.subplots(1, 4, figsize=(20, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
sns.histplot(adata.layers["log1p_norm_1e4"].sum(1), bins=100, kde=False, ax=axes[2])
axes[2].set_title("Shifted logarithm (defualt)")
sns.histplot(adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[3])
axes[3].set_title("log1p with Scran estimated size factors")
plt.show()

# Explanation:
# - Shifted Logarithm: Applies a log1p transformation to stabilize variance in counts.
# - Scran Normalization: Uses a clustering-based method to estimate size factors, 
#   which adjusts for differences in cell-specific biases before applying log1p normalization.

# In most cases, Shifted Logarithm Normalization is sufficient.

adata.write("./data/s4d8_quality_control_normalization.h5ad")
