
# -*- coding: utf-8 -*-

rm(list = ls())

install.packages("Seurat")

setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# BPCells
# Developed by Greenleaf Lab, BPCells is an R package that allows for computationally efficient single-cell analysis.
# It utilizes bit-packing compression to store counts matrices on disk and C++ code to cache operations.
# BPCells is an R package that allows for computationally efficient single-cell analysis.
# It utilizes bit-packing compression to store counts matrices on disk and C++ code to cache operations.

# presto
# Developed by Korunsky/Raychaudhari labs, presto performs a fast Wilcoxon rank sum test and auROC analysis.
# Seurat uses the presto package to perform fast differential expression.

# glmGamPoi
# Fit Gamma-Poisson Generalized Linear Models Reliably.

# install.packages("Signac")
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

# Signac
# Signac is designed for the analysis of single-cell chromatin data, including scATAC-seq.

# seurat-data
# SeuratData is a mechanism for distributing datasets in the form of Seurat objects using R's internal package and data management systems.

# azimuth
# Azimuth is a web application that uses an annotated reference dataset to automate the processing, analysis, and interpretation of a new single-cell RNA-seq experiment.

# seurat-wrappers
# SeuratWrappers is a collection of community-provided methods and extensions for Seurat, curated by the Satija Lab at NYGC. These methods comprise functionality not presently found in Seurat, and are able to be updated much more frequently.

library(Seurat)
library(SeuratData)
library(BPCells)

data_info <- SeuratData::AvailableData()
utils::View(data_info)

# PBMC 관련 데이터셋 (말초혈 단핵세포, Peripheral Blood Mononuclear Cells)

# Azimuth Reference


SeuratData::InstalledData()

SeuratData::InstallData("pbmc3k")
pbmc3k <- SeuratData::LoadData("pbmc3k")

pbmc3k
help(pbmc3k)

SeuratData::RemoveData("pbmc3k")
