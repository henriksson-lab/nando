% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{CellShapsToSeurat}
\alias{CellShapsToSeurat}
\title{Turn per-cell SHAPs into a seurat assay}
\usage{
CellShapsToSeurat(nando_dir)
}
\value{
Seurat assay, to be added by adata[\link{"SHAP"}] <- CellShapsToSeurat(nando_dir)
}
\description{
Turn per-cell SHAPs into a seurat assay
}
