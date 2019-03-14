## ---- echo=FALSE---------------------------------------------------------
#options(vignetteDocumentFormat=rmarkdown::all_output_formats("rwaveletvignette.Rmd"))

## ----Setup, include=FALSE------------------------------------------------
library(rwavelet)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5,
  fig.height=5, 
  fig.align="center"
)
set.seed(0)
par(mgp = c(2,0.5,0),
    mar = c(4,3,2,1),
    oma = c(1,1,0,0))

## ---- echo=TRUE----------------------------------------------------------
J <- 10; n <- 2^J; t <- (1:n) / n
name <- c('Bumps')
f <- MakeSignal(name, n)
plot(t, f, xlab="t", ylab="f(t)",
     type='l', lwd=1.2, main=name)

## ---- echo=TRUE----------------------------------------------------------
SNR <- 4
set.seed(1)
ssig <- sd(f)
sigma <- ssig / SNR
y <- f + rnorm(n, mean=0, sd=sigma)
plot(t, y, xlab="t", ylab="y", main=paste("Noisy", name))

## ---- echo=TRUE----------------------------------------------------------
qmf <- MakeONFilter('Daubechies', 10)
j0 <- 3
wc <- FWT_PO(y, j0, qmf)
wch <- wcs <- wcbs <- wcbh <- wc
wcf <- FWT_PO(f, j0, qmf)

## ---- echo=TRUE----------------------------------------------------------
PlotWaveCoeff(wc, j0, 0.5)
title("Noisy Wavelet Coefs")
PlotWaveCoeff(wcf, j0, 0.5)
title("Original Wavelet Coefs")

## ---- echo=TRUE----------------------------------------------------------
# estimate sigma using the Median Absolute Deviation
# using only the fine scale of wc
hatsigma <- MAD(wc[(2^(J-1)+1):2^J])
thr <- sqrt(2*log(length(y)))*hatsigma
# apply hard thresholding 
wch[(2^(j0)+1):n] <- HardThresh(wc[(2^(j0)+1):n], thr)
# plot the resulting coeficients estimates
PlotWaveCoeff(wch, j0, 0.5)
title("Estimated Wavelet Coefs")
fest <- IWT_PO(wch, j0, qmf)
snrout <- SNR(f, fest)
plot(t, fest, type='l', lwd=1.4, col='red', xlab="t", ylab="hat_f(t)",
	   main=format(round(snrout,2), nsmall=2))
matlines(t, f, type='l', lty=2)

## ------------------------------------------------------------------------
L <- 2 # block size
wcbh <- BlockThresh(wc, j0, hatsigma, L, qmf, thresh = "hard")
fest_bh <- IWT_PO(wcbh,j0,qmf)
PlotWaveCoeff(wcbh, j0, 0.5)
title("Estimated Wavelet Coefs")
plot(t, fest_bh, type='l', lwd=1.4, col='red', xlab="t", ylab="hat_f(t)",
	   main=format(round(SNR(f, fest_bh),2), nsmall=2))
matlines(t, f, type='l', lty=2)

## ---- echo=TRUE----------------------------------------------------------
wcs[(2^(j0)+1):n] <- SoftThresh(wc[(2^(j0)+1):n], thr)
PlotWaveCoeff(wcs, j0, 0.5)
title("Estimated Wavelet Coefs (soft)")
f_soft <- IWT_PO(wcs, j0, qmf)
plot(t, f_soft, type='l', lwd=1.4, col='red', xlab="t", ylab="hat_f(t)",
	   main=format(round(SNR(f, f_soft),2), nsmall=2))
matlines(t, f, type='l', lty=2)

## ------------------------------------------------------------------------
wcbs <- BlockThresh(wc, j0, hatsigma, L, qmf, thresh = "soft")
PlotWaveCoeff(wcbs, j0, 0.5)
fest_bs <- IWT_PO(wcbs,j0,qmf)
plot(t, fest_bs, type='l', lwd=1.4, col='red', xlab="t", ylab="hat_f(t)",
	   main=format(round(SNR(f, fest_bs),2), nsmall=2))
matlines(t, f, type='l', lty=2)

## ---- echo=TRUE, message=FALSE-------------------------------------------
library(imager)
name <- '../inst/extdata/lena.png'
f <- load.image(name)
plot(f, axes=F, interpolate=F, xlab="", ylab="")

## ---- echo=TRUE, message=FALSE-------------------------------------------
ssig <- sd(f)
sdnoise <- ssig/SNR
y <- f + rnorm(ncol(f)*nrow(f), mean=0, sd=sdnoise)
snrin <- SNR(f,y)
plot(y, axes=F, interpolate=F, xlab="", ylab="",
     main=format(round(snrin,2), nsmall = 2))

## ---- echo=TRUE, message=FALSE-------------------------------------------
wc <- FWT2_PO(as.array(squeeze(y)), j0, qmf)
wcb <- wc
thr <- 3*sdnoise
aT <- wc*(abs(wc)>thr)
fest <- IWT2_PO(aT, j0, qmf)
snrout <- SNR(f, fest)
plot(as.cimg(fest), axes=FALSE, xlab="", ylab="",
     main=format(round(snrout,2), nsmall=2))

## ---- eval=FALSE---------------------------------------------------------
#  library(misc3d)
#  library(rgl)
#  data("SLphantom")
#  n <- dim(SLphantom)[1]
#  contour3d(SLphantom,0)
#  rglwidget()

## ---- eval=FALSE---------------------------------------------------------
#  sig <- sd(SLphantom)
#  sdnoise <- sig/SNR
#  y <- SLphantom + rnorm(n^3, mean=0, sd=sdnoise)
#  wcf <- FWT3_PO(SLphantom,j0,qmf)
#  wc <- FWT3_PO(y,j0,qmf)
#  wcn <- wc
#  thr <- 3*sdnoise
#  wc <- wc*(abs(wc)>thr)
#  fhat <- IWT3_PO(wc, j0, qmf)

## ---- eval=FALSE---------------------------------------------------------
#  op <- par(mfrow=c(3,3), mgp = c(1.2, 0.5, 0), tcl = -0.2,
#  mar = .1 + c(0.1,0.1,0.1,0.1), oma = c(0,0,0,0))
#  ll=screen=list(z = 130, x = -80)
#  contour3d(SLphantom,0,color="gray", engine = "standard")
#  contour3d(y,0.05,color="gray", engine = "standard")
#  contour3d(fhat,0.2,color="gray", engine = "standard")
#  contour3d(wcf,0.1, color="gray", engine = "standard")
#  contour3d(wcn,0.05,color="gray", engine = "standard")
#  contour3d(wc,0.1,color="gray", engine = "standard")
#  plot(as.cimg(SLphantom[n/2,,]), axes=FALSE, xlab="", ylab="")
#  plot(as.cimg(y[n/2,,]), axes=FALSE, xlab="", ylab="")
#  plot(as.cimg(fhat[n/2,,]), axes=FALSE, xlab="", ylab="")

