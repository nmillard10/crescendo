# Crescendo
Crescendo for single-gene correction in single-cell RNA-sequencing data

Implementation of the Crescendo algorithm, which allows investigators to remove the effects of confounding factors by directly correcting gene expression count data.

## Installation

Install the current version of Crescendo from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/crescendo")
```

# Usage/Demos

## Quick-start code

The following code uses Crescendo to correct gene expression from a 10X dataset containing 3 batches (3pV1, 3pV2, and 5p)

``` r
library(crescendo)

# Load dataset with metadata, raw gene counts, and Harmony clusters (result from running the Harmony algorithm)
obj <- readRDS(system.file("extdata", "pbmc_obj.rds", package = "crescendo"))

# Set which genes to correct and parameters for coorrection (parameters explained in GettingStarted vignette)
genes_use <- c('TRAC', 'MS4A1')
prop <- 0.05
min_cells <- 50

batch_col <- 'batch'
id_col = 'CellID'
constant_umi <- TRUE
merge_clusters <- TRUE

mc.cores <- NULL
# mc.cores <- 1
alpha <- 0
lambda <- NULL

# Run Crescendo
corr_counts <- crescendo(
    Ycounts = obj$exprs_raw,
    meta_data = obj$meta_data,
    R = obj$R,
    genes_use = genes_use,
    prop = prop,
    min_cells = min_cells,
    batch_col = batch_col,
    id_col = id_col,
    constant_umi = TRUE,
    merge_clusters = TRUE,
    alpha = 0,
    lambda = NULL,
    seed = 2,
    return_obj = FALSE,
    mc.cores = mc.cores
)
corr_counts %>% str(1)
```
