% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsampling.r
\name{get_downsampled_idx}
\alias{get_downsampled_idx}
\title{Get the indices from downsampled data}
\usage{
get_downsampled_idx(
  meta_data,
  R,
  batch_var = "batch",
  prop = 0.1,
  min_cells = 50,
  return_idx = TRUE,
  clus_col = "sampled_clus",
  seed = 2,
  mc.cores = NULL,
  verbose = FALSE
)
}
\arguments{
\item{meta_data}{Metadata table that contains cell information (cell x info).}

\item{R}{A matrix of harmony clusters (cluster x cell) that are received as output from the Harmony algorithm.}

\item{batch_var}{Name of column containing batch information (e.g. batch, dataset, technology).}

\item{prop}{Fraction of cells to downsample to (e.g. 0.5 means downsample to arouond 50\% of original cells, 0.25 means 25\%, etc.).}

\item{min_cells}{Number of minimum cells to keep from each cell type in each batch.}

\item{return_idx}{Return just the indices of the downsampled metadata, or return the metadata table.}

\item{clus_col}{Name of column containing cluster information.}

\item{seed}{Seed for RNG.}

\item{mc.cores}{Number of cores to use for parallel functions.}

\item{verbose}{Print out progress messages.}
}
\value{
Indices of downsample data (or the whole metadata table if return_idx = FALSE)
}
\description{
To downsample data, we ensure that each cell-type in each batch is represented by a minimum number of cells (or if lower than
the minimum, all of the cells). This function can return just the indices, or the whole medatadata table.
}
