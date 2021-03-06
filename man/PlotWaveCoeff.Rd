% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotWaveCoeff.R
\name{PlotWaveCoeff}
\alias{PlotWaveCoeff}
\title{Spike-plot display of wavelet coefficients}
\usage{
PlotWaveCoeff(wc, L, scal)
}
\arguments{
\item{wc}{1-d wavelet transform.}

\item{L}{level of coarsest scale.}

\item{scal}{scale factor (0 ==> autoscale).}
}
\value{
A display of wavelet coefficients (coarsest level
        NOT included) by level and position.
}
\description{
Spike-plot display of wavelet coefficients
}
\examples{
x <- MakeSignal('Ramp', 128)
qmf <- MakeONFilter('Daubechies', 10)
L <- 3
scal <- 1
wc <- FWT_PO(x, L, qmf)
PlotWaveCoeff(wc,L,scal)
}
\seealso{
\code{\link{FWT_PO}}, \code{\link{IWT_PO}}, \code{\link{PlotSpikes}}.
}
