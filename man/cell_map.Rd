% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{cell_map}
\alias{cell_map}
\title{Plot heatmap of clustering coefficient matrix}
\usage{
cell_map(object, rank, main = "Cells", ...)
}
\arguments{
\item{object}{Object of class \code{scNMFSet}.}

\item{rank}{Rank value for which the cell map is to be displayed. 
The object must contain the corresponding slot: one element of 
          \code{coeff(object)[[k]]} for which 
          \code{ranks(object)[[k]]==rank}.}

\item{main}{Title of plot.}

\item{...}{Other arguments to be passed to \code{\link{heatmap}}, 
\code{\link{image}},
           and \code{\link{plot}}.}
}
\value{
\code{NULL}
}
\description{
Retrieve a coefficient matrix \code{H}
derived from factorization by rank value
and generate heatmap of its elements.
}
\examples{
set.seed(1)
x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
rownames(x) <- seq_len(10)
colnames(x) <- seq_len(100)
s <- scNMFSet(count=x,rowData=seq_len(10),colData=seq_len(100))
s <- vb_factorize(s,ranks=seq(2,5))
plot(s)
cell_map(s, rank=3)
}
