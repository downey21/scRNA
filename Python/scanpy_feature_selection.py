
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/feature_selection.html

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt

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
    filename="./data/s4d8_quality_control_normalization.h5ad",
)

ro.globalenv["adata"] = adata

# Method 1: deviance feature selection

# Higher deviance:
# The gene exhibits a structured (non-random) expression pattern with potential biological significance.
# It likely shows meaningful differences across cell populations.
# It is less likely to be noise and more likely to contain important biological information.

# Lower deviance:
# The gene expression is relatively stable and predictable.
# It shows low variability and may not differentiate between cell populations.
# It could be a candidate for removal in downstream analysis.

ro.r(
    """
    library(scry)
    library(SingleCellExperiment)
    sce <- scry::devianceFeatureSelection(adata, assay="X")
    binomial_deviance <- SingleCellExperiment::rowData(sce)$binomial_deviance
    """
)

# Retrieve binomial deviance from R to Python
binomial_deviance = np.array(ro.r("binomial_deviance")).T

# Select top 4000 features based on deviance
idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

# Store highly deviant features in adata.var
adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance

# Method 2: Compute highly variable genes using Scran normalization layer
sc.pp.highly_variable_genes(adata, layer="scran_normalization")

# Plot mean vs dispersion, highlighting highly deviant genes
ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()

adata.write("./data/s4d8_quality_control_normalization_feature_selection.h5ad")
