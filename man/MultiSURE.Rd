% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiSURE.R
\name{MultiSURE}
\alias{MultiSURE}
\title{Apply Shrinkage to Wavelet Coefficients}
\usage{
MultiSURE(wc, L)
}
\arguments{
\item{wc}{Wavelet Transform of noisy sequence with N(0,1) noise.}

\item{L}{low-frequency cutoff for Wavelet Transform.}
}
\value{
\code{ws} result of applying SUREThresh to each dyadic block.
}
\description{
SURE referes to Stein's Unbiased Risk Estimate.
}
