% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{meta_genes}
\alias{meta_genes}
\title{Find metagenes from basis matrix}
\usage{
meta_genes(object, rank, basis.matrix = NULL, max.per.cluster = 10,
  gene_names = NULL, subtract.mean = TRUE, log = TRUE)
}
\arguments{
\item{object}{Object of class \code{scNMFSet}.}

\item{rank}{Rank value for which metagenes are to be found.}

\item{basis.matrix}{Instead of an object containing basis 
matrices, the matrix itself can be provided.}

\item{max.per.cluster}{Maximum number of metagenes per cluster.}

\item{gene_names}{Names of genes to replace row names of basis matrix.}

\item{subtract.mean}{Standardize the matrix elements with means 
within each row.}

\item{log}{Use geometric mean and division instead of arithmetic 
mean and subtraction with \code{subtract.mean}.}
}
\value{
List of vectors each containing metagene names of clusters.
}
\description{
Retrieve a basis matrix from an object and find metagenes.
}
\examples{
set.seed(1)
x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
rownames(x) <- seq_len(10)
colnames(x) <- seq_len(100)
s <- scNMFSet(count=x,rowData=seq_len(10),colData=seq_len(100))
s <- vb_factorize(s,ranks=seq(2,5))
meta_genes(s, rank=4)
}
