% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{create_design_matrix}
\alias{create_design_matrix}
\title{Create a design matrix.}
\usage{
create_design_matrix(meta_data, R)
}
\arguments{
\item{meta_data}{Metadata table that contains cell information (cell x info).}

\item{R}{A matrix containing cluster assignments for each cell (cluster x cell).}
}
\value{
A design matrix that contains batch and cluster membership information.
}
\description{
This creates a desgin matrix from a meta-data table and cluster identities.
}
