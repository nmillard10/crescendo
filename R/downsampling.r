#' Sample a discrete harmony cluster based off of soft-cluster probabilities.
#'
#' Harmony soft-clusters give probabilities for each cell. In order to downsample data, we sample a discrete cluster for each cell
#' based on its harmony cluster probabilities.
#'
#' @param harmony_clusters A matrix of harmony clusters (cluster x cell) that are received as output from the Harmony algorithm.
#' @param seed Seed for RNG.
#' @param mc.cores Number of cores to use for sampling (useful for large datasets).
#' @return A vector containing a discrete cluster assignment for each cell design matrix that contains batch and cluster membership 
#' information.
#' @importFrom tidyselect all_of
#' @export
sample_harmony_clus <- function(harmony_clusters, seed = NULL, mc.cores = NULL){
    RNGkind("L'Ecuyer-CMRG")
    if(is.null(mc.cores)){
        mc.cores <- max(1, detectCores() - 1)
    }
    if(!is.null(seed)){
        set.seed(seed)
    }
    sampled_clus <- factor(mclapply(seq(nrow(harmony_clusters)), function(x){
        sample(x = colnames(harmony_clusters), size = 1, prob = harmony_clusters[x, ])
    }, mc.cores = mc.cores, mc.preschedule = TRUE), levels = colnames(harmony_clusters))
    return(sampled_clus)
}

#' Downsample data in a batch and cell-type aware manner.
#'
#' To downsample data, we ensure that each cell-type in each batch is represented by a minimum number of cells (or if lower than
#' the minimum, all of the cells).
#'
#' @param meta_data Metadata table that contains cell information (cell x info).
#' @param id_col Name of column containing unique cell IDs.
#' @param batch_col Name of column containing batch information (e.g. batch, dataset, technology).
#' @param clus_col Name of column containing cluster information.
#' @param prop Fraction of cells to downsample to (e.g. 0.5 means downsample to arouond 50\% of original cells, 0.25 means 25\%, etc.).
#' @param min_cells Number of minimum cells to keep from each cell type in each batch.
#' @param seed Seed for RNG.
#' @return A metadata table containing the downsampled cells.
#' @importFrom dplyr slice_sample
#' @export
downsample_batch_clus <- function(
    meta_data, 
    id_col = 'CellID', 
    batch_col = 'batch', 
    clus_col = 'sampled_clus', 
    prop,
    min_cells,
    seed = NULL
){
    col_get <- c(id_col, batch_col, clus_col)
    temp_meta <- meta_data %>% dplyr::select(all_of(col_get))
    temp_meta$idx <- seq(nrow(temp_meta))
    temp_meta <- split(temp_meta, temp_meta[[batch_col]])
    temp_meta <- lapply(temp_meta, function(x){
        split(x, x[[clus_col]])
    })
    
    if(!is.null(seed)){
        set.seed(seed)
    }
    temp_meta <- lapply(temp_meta, function(x){
        lapply(x, function(y){
            if(nrow(y) > min_cells){
                samp_n <- max(min_cells, ceiling(nrow(y) * prop))
                new_frame <- y %>% dplyr::slice_sample(n = samp_n)
                return(new_frame)
            } else{
                return(y)
            }
        })
    })
    
    temp_meta <- unlist(temp_meta, recursive = FALSE)
    temp_meta <- rbindlist(temp_meta)
    return(temp_meta)
}