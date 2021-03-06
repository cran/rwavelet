% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block_partition2d.R
\name{block_partition2d}
\alias{block_partition2d}
\title{Construct 2d block partition}
\usage{
block_partition2d(x, L)
}
\arguments{
\item{x}{(noisy) wc at a given scale.}

\item{L}{block size.}
}
\value{
\code{out} partition of coefficients by block.
}
\description{
Group the coefficients into blocks (or groups) of size L.
}
\examples{
x <- matrix(rnorm(2^2), ncol=2)
j0 <- 0
qmf <- MakeONFilter('Haar')
wc <- FWT2_PO(x, j0, qmf)
L <- 2
wcb <- block_partition2d(wc, L)
}
\seealso{
\code{\link{invblock_partition2d}}
}
