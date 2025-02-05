
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0 # 0, 1, 2 (default), 3
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

#  single-cell multiomics data from bone marrow mononuclear cells of
# 12 healthy human donors measured at four different sites to obtain nested batch effects
# In this tutorial, we will use one batch of the aforementioned dataset, sample 4 of donor 8,
# to showcase the best practices for scRNA-seq data preprocessing.

# Read 10x Genomics HDF5 file
# Output: AnnData object
adata = sc.read_10x_h5(
    filename="./data/filtered_feature_bc_matrix.h5",
    backup_url="https://figshare.com/ndownloader/files/39546196",
)
adata # n_obs: number of cells; n_vars: number of genes; var: meta

adata.var_names_make_unique()
adata

adata.obs.head() # meta data of cells
adata.var.head() # meta data of genes

adata.X[:5, :5]
adata.X[:5, :5].todense()
adata.X[:5, :5].toarray()
pd.DataFrame(
    adata[:5, :5].X.toarray(),
    index=adata.obs_names[:5],
    columns=adata.var_names[:5]
)

adata[:, 0].X
adata[:, 0].X.toarray()

adata[0, :].X
adata[0, :].X.toarray()

# Filtering low quality cells

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
adata.var.head()

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
adata.obs.head()
adata.var.head()

p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
p1.figure.show()
sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.head()
adata.obs.outlier.value_counts()

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.head()
adata.obs.mt_outlier.value_counts()

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Correction of ambient RNA
