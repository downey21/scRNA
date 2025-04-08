# Single Cell Data Preprocessing

## Raw (Read) Data

### Raw Data Generation and Processing
- Sequencing machines generate raw data that must be processed into a structured format.
- Data can be represented as count matrices when unique molecular identifiers (UMIs) are used.

### Raw Data Processing Pipelines
Several pipelines are available for processing raw sequencing data:
- **Cell Ranger**
- **indrops**
- **SEQC**
- **zUMIs**

### Key Processing Steps
1. Quality Control (QC) of read data
2. Demultiplexing: Assigning reads to specific barcodes linked to individual cells.
3. Genome alignment: Mapping reads to a reference genome.
4. Quantification: Counting the number of transcripts for each gene in each cell.

### Understanding Count Matrices and Barcodes
- The output of processing is typically a count matrix with dimensions: Number of barcodes × Number of transcripts.
- Why use barcodes instead of cells?
  - In single-cell RNA sequencing, barcodes are assigned to individual mRNA molecules rather than directly identifying cells.
  - Not all barcodes correspond to a single cell due to technical limitations.
- Possible errors in barcode assignment:
  - Doublets: A single barcode mistakenly tags multiple cells, causing mixed gene expression profiles.
  - Empty droplets/wells: A barcode is assigned to an empty droplet or well, meaning no actual cell is present.

## Cell-level QC (Filtering Cells)

### Purpose
- Ensures that only viable cells are used for downstream analysis.
- Removes low-quality, dying, or damaged cells.

### Key QC Metrics
Three primary QC covariates are used to assess cell quality:
- Count depth: Total number of counts per barcode.
- Number of genes detected: Helps identify low-quality cells.
- Mitochondrial fraction: High levels indicate cell damage or apoptosis.

### Outlier (Cell) Detection: Filtering Cells
Best practices: Use multiple QC metrics together rather than filtering based on a single threshold.

- Low count depth, few genes, high mitochondrial fraction: Dying or damaged cells.
- High count depth, many genes: Doublets (barcodes mistakenly tagging multiple cells).
- Tools for detecting doublets: DoubletDecon, Scrublet, Doublet Finder.

Tips:
- Cells with high mitochondrial reads may be biologically active, not just damaged.
- High-count cells could be larger cells, not necessarily doublets.

### Handling Heterogeneous Datasets
- Different cell types may have different QC distributions.
- Set adaptive thresholds based on dataset characteristics.
- Avoid filtering too aggressively: Some low-quality-looking cells might be valid subpopulations.

## Gene-level QC (Filtering Genes)

### Filtering Genes
- Remove genes expressed in very few cells.
- Be cautious in high-dropout datasets: Strict filtering may remove real biological signals.

### Correcting Ambient RNA Contamination
- Ambient RNA comes from lysed cells, distorting analysis.
- Use **SoupX** or other correction methods to filter out contaminating signals.

## Iterative QC Optimization
- No fixed QC thresholds work for all datasets; adjustments may be needed.
- Start with permissive QC thresholds, then refine based on downstream analysis results.
- Avoid modifying thresholds just to improve statistical test results (data peeking).
- Use visualizations to validate QC decisions before finalizing filters.

## Normalization

### Purpose
- Normalization adjusts for differences in count depth caused by sequencing variability.
- Ensures that statistical methods can be more accurately applied in subsequent analyses.

### Normalization Methods
(No single method is universally optimal)
- **Shifted Logarithm Normalization**
- **Scran’s Pooling-Based Size Factor Estimation Method and log(x+1) transformation**

## Batch Effects Correction

### Batch Effects
- Occur when cells are processed under different experimental conditions (e.g., different sequencing lanes, time points).
- Can distort gene expression comparisons and mask biological signals.

### Batch Correction Methods
- **ComBat**
- Preventing Batch Effects is Ideal (**Good experimental design**)

## Expression Recovery (Denoising/Imputation)

- Dropouts (zero counts due to sampling noise) are common in scRNA-seq data.
- Expression recovery estimates missing values to improve gene–gene correlation analysis.

### Common Imputation Tools
- **MAGIC**
- **DCA**
- **scVI**
- **SAVER**
- **scImpute**

### Potential Issues
- No consensus on whether imputation should be used.
- Scalability remains a challenge for large datasets.

## Feature Selection

### Purpose
- Reduces dataset dimensionality by keeping only informative genes.
- Highly variable genes (HVGs) are commonly used for analysis.

### Best Practice
- Perform feature selection after technical data correction to avoid selecting genes influenced by other effects.

## Dimensionality Reduction

### Purpose
- Further reduces the dimensionality of single-cell expression matrices after feature selection.
- Embeds the dataset into a low-dimensional space that captures the underlying biological structure.
- Maps the dataset into 2D or 3D for visualization.

### Methods
- **PCA** (Principal Component Analysis)
- **t-SNE** (t-distributed Stochastic Neighbor Embedding)
- **UMAP** (Uniform Manifold Approximation and Projection)

## References
This README is based on the followings:

- Heumos, Lukas, et al. "Best practices for single-cell analysis across modalities." *Nature Reviews Genetics* 24.8 (2023): 550-572.

- Luecken, Malte D., and Fabian J. Theis. "Current best practices in single-cell RNA-seq analysis: a tutorial." *Molecular systems biology* 15.6 (2019): e8746.

- [The Single-Cell Best Practices online book](https://www.sc-best-practices.org/preamble.html).
