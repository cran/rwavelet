% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BlockThresh.R
\name{BlockThresh}
\alias{BlockThresh}
\title{1d wavelet Block Thresholding}
\usage{
BlockThresh(wc, j0, hatsigma, L, qmf, thresh = "hard")
}
\arguments{
\item{wc}{wavelet coefficients.}

\item{j0}{coarsest decomposition scale.}

\item{hatsigma}{estimator of noise variance.}

\item{L}{Block size (n mod L must be 0).}

\item{qmf}{Orthonormal quadrature mirror filter.}

\item{thresh}{'hard' or 'soft'.}
}
\value{
\code{wcb} wavelet coefficient estimators.
}
\description{
This function is used for thresholding coefficients by group (or block)
according to the hard or soft thresholding rule.
}
\examples{
n <- 64
x <- MakeSignal('Ramp', n)
sig <- 0.01
y <- x + rnorm(n, sd=sig)
j0 <- 1
qmf <- MakeONFilter('Daubechies',8)
wc <- FWT_PO(y, j0, qmf)
L <- 2
wcb <- BlockThresh(wc, j0, sig, L, qmf, "hard")
}
\seealso{
\code{\link{invblock_partition}}, \code{\link{invblock_partition}}.
}
