
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

# Read 10x Genomics HDF5 file (.h5)
# Output: AnnData object (Python Scanpy .h5ad)
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

# The number of counts per barcode (count depth)
# The number of genes per barcode
# The fraction of top 20 genes
# The fraction of counts from mitochondrial genes per barcode

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

# p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
# p1.figure.show()
# sc.pl.violin(adata, 'total_counts')
# p2 = sc.pl.violin(adata, "pct_counts_mt")
# p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# MAD = median(|X_i - median(X)|)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5) # The number of counts per barcode (count depth)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5) # The number of genes per barcode
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5) # The fraction of top 20 genes
)
adata.obs.head()
adata.obs.outlier.value_counts()

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8 # The fraction of counts from mitochondrial genes per barcode
)
adata.obs.head()
adata.obs.mt_outlier.value_counts()

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

# sc.pl.violin(adata, 'total_counts')
# sc.pl.violin(adata, "pct_counts_mt")
# p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Correction of ambient RNA

# Cell-free mRNA molecules, also known as ambient RNA,
# can confound the number of observed counts and can be seen as background contamination.
# It is important to correct droplet-based scRNA-seq datasets for cell-free mRNA
# as it may distort the interpretation of the data in our downstream analysis.

# SoupX does not directly filter genes but adjusts their expression values.
# However, after applying sc.pp.filter_genes(adata, min_cells=20), lowly expressed genes are removed.

import anndata2ri # AnnData (.h5ad) -> SingleCellExperiment (.RData)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

pandas2ri.activate()
anndata2ri.activate()

adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp) # set total count equal (e.g., 10000) for all cells
sc.pp.log1p(adata_pp)

# Cell Clustering for SoupX
sc.pp.pca(adata_pp) # PCA for Leiden Clustering
print(adata_pp.obsm["X_pca"].shape) # (n_cells, n_pcs)
print(adata_pp.uns["pca"]) # PCA meta
print(adata_pp.varm["PCs"].shape)  # (n_genes, n_pcs)

sc.pp.neighbors(adata_pp) # Nearest Neighbors for Leiden Clustering
print(adata_pp.obsp["distances"]) # distance matrix
print(adata_pp.obsp["connectivities"]) # adjacency matrix
print(adata_pp.uns["neighbors"])  # parameters

# sc.tl.leiden(adata_pp, key_added="soupx_groups") # Leiden Clustering
sc.tl.leiden(adata_pp, key_added="soupx_groups", flavor="igraph", n_iterations=2, directed=False) # Leiden Clustering
print(adata_pp.obs["soupx_groups"].value_counts())

# Preprocess variables for SoupX
soupx_groups = adata_pp.obs["soupx_groups"]

del adata_pp

cells = adata.obs_names
genes = adata.var_names
data = adata.X.T

adata_raw = sc.read_10x_h5(
    filename="./data/raw_feature_bc_matrix.h5",
    backup_url="https://figshare.com/ndownloader/files/39546217",
)
adata_raw.var_names_make_unique()
data_tod = adata_raw.X.T

del adata_raw

ro.r('library(SoupX)')

ro.globalenv['data'] = data
ro.globalenv['data_tod'] = data_tod
ro.globalenv['genes'] = genes
ro.globalenv['cells'] = cells
ro.globalenv['soupx_groups'] = soupx_groups

ro.r('''
# specify row and column names of data
rownames(data) <- genes
colnames(data) <- cells

# ensure correct sparse format for table of counts and table of droplets
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Generate SoupChannel Object for SoupX 
sc <- SoupX::SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf <- data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc <- SoupX::setSoupProfile(sc, soupProf)

# Set cluster information in SoupChannel
sc <- SoupX::setClusters(sc, soupx_groups)

# Estimate contamination fraction
sc <- SoupX::autoEstCont(sc, doPlot=FALSE)

# Infer corrected table of counts and rount to integer
out <- SoupX::adjustCounts(sc, roundToInt = TRUE)
''')

out = ro.r('out')

adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
adata.X = adata.layers["soupX_counts"]

print(f"Total number of genes: {adata.n_vars}")

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")

# Doublet Detection

# In scRNA-seq (single-cell RNA sequencing), each cell should have a unique barcode.
# However, in some cases, multiple cells may be captured under the same barcode and sequenced together.
# This results in RNA from two cells being mixed, making them appear as a single cell.
# This phenomenon is called a doublet, and tools like scDblFinder are used to detect it.

# When multiple batches are combined, batch effects can cause normal cells to be mistakenly identified as doublets.
# Therefore, it is recommended to perform doublet detection separately for each batch before merging the data.

ro.r('''
library(Seurat)
library(scater)
library(scDblFinder) # R package for doublet detection
library(BiocParallel)
''')

data_mat = adata.X.T
ro.globalenv["data_mat"] = data_mat

ro.r('''
set.seed(123)
sce <- scDblFinder(SingleCellExperiment(list(counts=data_mat)))
doublet_score <- sce$scDblFinder.score
doublet_class <- sce$scDblFinder.class
''')

doublet_score = ro.r("doublet_score")
doublet_class = ro.r("doublet_class")

adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class
adata.obs.scDblFinder_class.value_counts()

adata.write("./data/s4d8_quality_control.h5ad")
# adata = sc.read_h5ad("./data/s4d8_quality_control.h5ad")
# library(zellkonverter); zellkonverter::writeH5AD(adata, "./data/s4d8_quality_control.h5ad"); adata <- zellkonverter::readH5AD("./data/s4d8_quality_control.h5ad")

# Seurat::SaveH5Seurat(seurat_object, filename = "single_cell_data.h5seurat")
# seurat_object <- Seurat::LoadH5Seurat("single_cell_data.h5seurat")

# ro.r('''
# library(SeuratDisk)

# # .h5seurat -> .h5ad
# SeuratDisk::Convert("single_cell_data.h5seurat", dest = "h5ad", overwrite = TRUE)

# # .h5ad -> .h5seurat
# SeuratDisk::Convert("single_cell_data.h5ad", dest = "h5seurat", overwrite = TRUE)
# ''')

print("Before filtering:")
print(f"Total number of cells before filtering: {adata.n_obs}")
print(adata.obs["scDblFinder_class"].value_counts())

adata = adata[adata.obs["scDblFinder_class"] == "singlet"].copy()

print("After filtering:")
print(f"Total number of cells after filtering: {adata.n_obs}")
print(adata.obs["scDblFinder_class"].value_counts())
