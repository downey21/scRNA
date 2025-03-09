
# -*- coding: utf-8 -*-

rm(list = ls())

library(Seurat)
library(SeuratData)
library(dplyr)
library(patchwork)

# pbmc.data <- Seurat::Read10X(data.dir = "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/")
# pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# SeuratData::InstallData("pbmc3k")
# pbmc <- SeuratData::LoadData("pbmc3k")

# pbmc
# # help(pbmc3k)

url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
curl::curl_download(url = url, destfile = file.path("./data", basename(url)))
untar(tarfile = file.path("./data", basename(url)), exdir = "./data")

url_raw <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz"
curl::curl_download(url = url_raw, destfile = file.path("./data", basename(url_raw)))
untar(tarfile = file.path("./data", basename(url_raw)), exdir = "./data")

pbmc_data <- Seurat::Read10X(data.dir = "./data/filtered_gene_bc_matrices/hg19/")
rownames(pbmc_data) <- gsub("_", "-", rownames(pbmc_data))
pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data, project = "pbmc3k")

pbmc_raw_data <- Seurat::Read10X(data.dir = "./data/raw_gene_bc_matrices/hg19/")
rownames(pbmc_raw_data) <- gsub("_", "-", rownames(pbmc_raw_data))
pbmc_raw <- Seurat::CreateSeuratObject(counts = pbmc_raw_data, project = "pbmc3k_raw")

head(pbmc@meta.data)

# table(pbmc$seurat_annotations)
# table(pbmc@meta.data$seurat_annotations)

# Filtering low quality cells

# The number of counts per barcode (count depth)
# The number of genes per barcode
# The fraction of top 20 genes
# The fraction of counts from mitochondrial genes per barcode

# The [[ ]] operator can add columns to object metadata.
# Compute Mitochondrial, Ribosomal, Hemoglobin Gene Percentage
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.ribo"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL")
pbmc[["percent.hb"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^HB[^(P)]")

head(pbmc@meta.data)

# Visualize QC metrics as a violin plot
# Seurat::VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# plot1 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# Outlier Removal using MAD
mad_outliers <- function(x, nmads = 5) {
    median_x <- median(x, na.rm = TRUE)
    mad_x <- mad(x, na.rm = TRUE)

    return(x < (median_x - nmads * mad_x) | x > (median_x + nmads * mad_x))
}

calculate_top_20_pct <- function(seurat_obj) {
    counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")

    top_20_counts <- apply(counts, 2, function(x) {
        sorted_x <- sort(x, decreasing = TRUE)
        top_20_sum <- sum(sorted_x[1:min(20, length(sorted_x))], na.rm = TRUE)

        return(top_20_sum / sum(x, na.rm = TRUE) * 100)
    })
    return(top_20_counts)
}

pbmc[["pct_counts_in_top_20_genes"]] <- calculate_top_20_pct(pbmc)

pbmc[["log1p_nFeature_RNA"]] <- log1p(pbmc[["nFeature_RNA"]])
pbmc[["log1p_nCount_RNA"]] <- log1p(pbmc[["nCount_RNA"]])

pbmc[["outlier"]] <- (
    mad_outliers(pbmc[["log1p_nFeature_RNA"]][,1]) |
    mad_outliers(pbmc[["log1p_nCount_RNA"]][,1]) |
    mad_outliers(pbmc[["pct_counts_in_top_20_genes"]][,1])
)

pbmc[["mt_outlier"]] <- (
    mad_outliers(pbmc[["percent.mt"]][,1], nmads = 3) | (pbmc[["percent.mt"]][,1] > 8)
)

dim(pbmc)
dim(pbmc@meta.data)
head(pbmc@meta.data)

pbmc <- pbmc[, !(pbmc@meta.data$outlier) & !(pbmc@meta.data$mt_outlier)]

dim(pbmc)
dim(pbmc@meta.data)

str(pbmc)
head(pbmc@assays$RNA@meta.data)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowMeans(GetAssayData(pbmc, layer = "counts")),
    col.name = "mean_counts"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowSums(GetAssayData(pbmc, layer = "counts") > 0),
    col.name = "n_cells"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowSums(GetAssayData(pbmc, layer = "counts")),
    col.name = "total_counts"
)

rownames(pbmc@assays$RNA@meta.data) <- rownames(pbmc)

head(pbmc@assays$RNA@meta.data)

# Ambient RNA Correction using SoupX

# Cell-free mRNA molecules, also known as ambient RNA,
# can confound the number of observed counts and can be seen as background contamination.
# It is important to correct droplet-based scRNA-seq datasets for cell-free mRNA
# as it may distort the interpretation of the data in our downstream analysis.

library(scater)
library(scran)
library(igraph)
library(SoupX)

sce <- Seurat::as.SingleCellExperiment(pbmc)

sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
snn_graph <- scran::buildSNNGraph(sce)
sce$clusters <- igraph::cluster_louvain(snn_graph)$membership

raw_counts <- Seurat::GetAssayData(pbmc_raw, assay = "RNA", layer = "counts")
filtered_counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")

soupX_obj <- SoupX::SoupChannel(raw_counts, filtered_counts)
soupX_obj <- SoupX::setClusters(soupX_obj, sce$clusters)
soupX_obj <- SoupX::autoEstCont(soupX_obj, doPlot = FALSE)

corrected_counts <- SoupX::adjustCounts(soupX_obj, roundToInt = TRUE)

pbmc[["RNA_original"]] <- Seurat::CreateAssayObject(counts = Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts"))
pbmc@meta.data$nCount_RNA_original <- NULL
pbmc@meta.data$nFeature_RNA_original <- NULL
pbmc <- Seurat::SetAssayData(pbmc, assay = "RNA", layer = "counts", new.data = corrected_counts)
SeuratObject::Key(pbmc[["RNA"]]) <- "rna_"
Seurat::DefaultAssay(pbmc) <- "RNA"

# Filter Lowly Expressed Genes
counts_matrix <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
genes_to_keep <- rowSums(counts_matrix > 0) >= 20
pbmc <- subset(pbmc, features = rownames(counts_matrix)[genes_to_keep])
pbmc@meta.data$nCount_RNA_original <- NULL
pbmc@meta.data$nFeature_RNA_original <- NULL

counts_matrix <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
nCount_RNA <- Matrix::colSums(counts_matrix)
nFeature_RNA <- Matrix::colSums(counts_matrix > 0)
pbmc@meta.data$nCount_RNA <- nCount_RNA[rownames(pbmc@meta.data)]
pbmc@meta.data$nFeature_RNA <- nFeature_RNA[rownames(pbmc@meta.data)]
pbmc@meta.data$log1p_nFeature_RNA <- log1p(pbmc@meta.data$nFeature_RNA)
pbmc@meta.data$log1p_nCount_RNA <- log1p(pbmc@meta.data$nCount_RNA)

# Doublet Detection using scDblFinder

# In scRNA-seq (single-cell RNA sequencing), each cell should have a unique barcode.
# However, in some cases, multiple cells may be captured under the same barcode and sequenced together.
# This results in RNA from two cells being mixed, making them appear as a single cell.
# This phenomenon is called a doublet, and tools like scDblFinder are used to detect it.

# When multiple batches are combined, batch effects can cause normal cells to be mistakenly identified as doublets.
# Therefore, it is recommended to perform doublet detection separately for each batch before merging the data.

library(scDblFinder)

sce <- Seurat::as.SingleCellExperiment(pbmc)
sce <- scDblFinder::scDblFinder(sce)

pbmc[["scDblFinder_score"]] <- sce$scDblFinder.score
pbmc[["scDblFinder_class"]] <- sce$scDblFinder.class

# Remove doublets
pbmc <- subset(pbmc, subset = scDblFinder_class == "singlet")

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowMeans(GetAssayData(pbmc, layer = "counts")),
    col.name = "mean_counts"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowSums(GetAssayData(pbmc, layer = "counts") > 0),
    col.name = "n_cells"
)

pbmc[["RNA"]] <- Seurat::AddMetaData(pbmc[["RNA"]],
    metadata = rowSums(GetAssayData(pbmc, layer = "counts")),
    col.name = "total_counts"
)

rownames(pbmc@assays$RNA@meta.data) <- rownames(pbmc)

base::saveRDS(pbmc, file = "./data/pbmc3k_qc.rds")
# pbmc <- base::readRDS(file = "./data/pbmc3k_qc.rds")
