
# -*- coding: utf-8 -*-

# https://www.sc-best-practices.org/cellular_structure/clustering.html

import scanpy as sc

sc.settings.verbosity = 0  # 0, 1, 2 (default), 3
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

adata = sc.read(
    filename="./data/s4d8_quality_control_normalization_feature_selection_dimensionality_reduction.h5ad"
)

sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

sc.pl.umap(
    adata,
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
)
