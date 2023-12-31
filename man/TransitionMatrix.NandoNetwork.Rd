% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{TransitionMatrix.NandoNetwork}
\alias{TransitionMatrix.NandoNetwork}
\title{Make a transition matrix}
\usage{
\method{TransitionMatrix}{NandoNetwork}(
  nandonet,
  subsetGene = NULL,
  subsetTF = NULL,
  includegenes = c()
)
}
\arguments{
\item{subsetGene, subsetTF}{Subset of genes to include. Default is all}

\item{includegenes}{Names of genes to forcefully add as row/columns in the matrix}

\item{onegrn}{A NandoNet}
}
\value{
Transition matrix, in sparse matrix format to save memory
}
\description{
Make a transition matrix
}
