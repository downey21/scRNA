
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/preprocessing_visualization/normalization.html

# Normalization

# It is the process of adjusting the variability caused by sampling effects to bring the data variance within a consistent range.
# This ensures that statistical methods can be more accurately applied in subsequent analyses.

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.sparse import issparse

import anndata2ri # AnnData (.h5ad) -> SingleCellExperiment (.RData)
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

sc.settings.verbosity = 0 # 0, 1, 2 (default), 3
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

pandas2ri.activate()
anndata2ri.activate()

adata = sc.read(
    filename="./data/download_s4d8_quality_control.h5ad",
    backup_url="https://figshare.com/ndownloader/files/40014331",
)

p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False)

# 
