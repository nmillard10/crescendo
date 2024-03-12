#' Implementation of the Crescendo algorithm.
#'
#' Runs the Crescendo algorithm by performing the downsampling (if necessary), estiimation, marginalization, and matching steps.
#'
#' @param Ycounts A matrix containing un-normalized gene count information (cells x genes).
#' @param meta_data Metadata table that contains cell information (cell x info).
#' @param R A matrix containing cluster assignments for each cell (cluster x cell). 
#' @param genes_use Genes to estimate. If none supplied, will automatically estimate all genes in Ycounts.
#' @param id_col Name of column containing unique cell IDs.
#' @param batch_col Name of column containing batch information (e.g. batch, dataset, technology).
#' @param prop Fraction of cells to downsample to (e.g. 0.5 means downsample to arouond 50\% of original cells, 0.25 means 25\%, etc.).
#' @param min_cells Number of minimum cells to keep from each cell type in each batch.
#' @param constant_umi Determines whether to sample counts at a constant UMI level vs. variable (recommended to be constant).
#' @param merge_clusters Determines whether to merge redundant Harmony clusters (merging decreases computational costs).
#' @param lambda Regularization parameter.
#' @param alpha Elastic net regularization parameter. 1 indicates lasso penalty, 0 indicates ridge penalty.
#' @param seed Seed for RNG.
#' @param return_obj Return only corrected counts, or results from each step.
#' @param mc.cores Number of cores to use for sampling (useful for large datasets).
#' @param impute_reads If not NULL, supply a vector of UMI to perform gene imputation in a higher UMI context than observed in the data.
#' @return A dataframe containing corrected gene count expressions, or a list containing results from each Crescendo step.
#' @import tidyverse
#' @import magrittr
#' @import parallel
#' @import data.table
#' @import Matrix
#' @import glmnet
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats coef median model.matrix sd
#' @export
crescendo <- function(Ycounts, meta_data, R, genes_use = NULL, id_col = 'CellID', batch_col = 'batch',
                      prop = 1, min_cells = 50, constant_umi = TRUE, merge_clusters = TRUE,
                      lambda = NULL, alpha = 0, seed = 2, return_obj = FALSE, mc.cores = NULL, impute_reads = NULL){
    if(!(batch_col %in% colnames(meta_data))){
        stop('Provided column does not exist. Please input name of batch column')
    }
    if(!(id_col %in% colnames(meta_data))){
        stop('Provided column does not exist. Please input name of cell ID column')
    }
    
    if(is.null(mc.cores)){
        mc.cores <- max(1, detectCores() - 1)
    }
    meta_data$batch <- meta_data[[batch_col]]
    
    # Merge redundant clusters
    if(merge_clusters){
        R <- t(merge_redundant_clusters(R, 0.99))
        colnames(R) <- paste0('R', 1:ncol(R))
    } else{
        R <- t(R)
        colnames(R) <- paste0('R', 1:ncol(R))
    }
    
    # Use selected genes is applicable
    if(!is.null(genes_use)){
        Ycounts <- Ycounts[genes_use, , drop = FALSE]
    } else{
        genes_use <- rownames(Ycounts)
    }
    Ycounts <- t(Ycounts)
    
    offset <- log(meta_data$nUMI)
    if(constant_umi){
        final_reads <- log(rep(median(meta_data$nUMI), length(meta_data$nUMI)))
    } else{
        final_reads <- offset
    }
    if(!is.null(impute_reads)){
        final_reads <- impute_reads
    }
    
    message('Creating design matrix')
    design_full <- create_design_matrix(meta_data = meta_data, R = R, id_col = id_col)
    # Designate terms to keep (e.g. cluster terms)
    terms_keep <- grep('^R\\d+$|disease|(Intercept)', colnames(design_full), value = TRUE)
    # Designate terms to remove (e.g. interaction terms)
    terms_remove_crossed <- terms_remove_crossed <- grep(':batch', colnames(design_full), value = TRUE)
    
    # Downsampling
    if(prop != 1){
        message('Getting discrete Harmony cluster')
        meta_data$sampled_clus <- sample_harmony_clus(
            harmony_clusters = R,
            seed = seed,
            mc.cores = mc.cores
        )
        message('Downsampling')
        temp_meta <- downsample_batch_clus(
            meta_data = meta_data,
            id_col = id_col,
            batch_col = batch_col,
            clus_col = 'sampled_clus',
            min_cells = min_cells,
            prop = prop,
            seed = seed
        )

        temp_Ycounts <- Ycounts[temp_meta$idx, , drop = FALSE]
        temp_R <- R[temp_meta$idx, ]
        temp_offset <- offset[temp_meta$idx]
        temp_design <- design_full[temp_meta$idx, , drop = FALSE]
    } else{
        message('No downsampling selected')
        temp_meta <- meta_data
        temp_Ycounts <- Ycounts
        temp_R <- R
        temp_offset <- offset
        temp_design <- design_full
    }
    
    # Estimation
    message('Estimating')
    betas <- estimate_betas(
        Ycounts = temp_Ycounts,
        R = temp_R,
        design_full = temp_design,
        lambda = lambda,
        alpha = alpha,
        genes_use = genes_use,
        offset = temp_offset,
        mc.cores = mc.cores
    )
    
    message('Marginalizing')
    expected_counts_marginalized <- marginalize_effects(
        R = R, 
        design_full = design_full, 
        betas = betas, 
        terms_keep = terms_keep, 
        terms_remove_crossed = terms_remove_crossed, 
        offset = final_reads
    )

    message('Matching')
    corrected_matching <- match_cdfs(
        expected_counts_marginalized = expected_counts_marginalized, 
        Ycounts = Ycounts, 
        design_full = design_full, 
        betas = betas, 
        offset = final_reads
    )

    if (return_obj) {
        list(
            betas = betas, ## from estimation
            expected_counts_marginalized = expected_counts_marginalized, ## from marginalization
            corrected_matching = corrected_matching ## from matching
        )
    } else {
        return(corrected_matching)            
    }
    
}