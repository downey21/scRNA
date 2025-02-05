
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

import numpy as np
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

# Filtering low quality cells
