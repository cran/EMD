\name{emd2d}
\alias{emd2d}
\title{Bidimenasional Empirical Mode Decomposition}
\description{
This function performs the bidimenasional empirical mode decomposition
utilizing extrema detection based
on the equivalence relation between neighboring pixels.
}
\usage{
emd2d(z, x = NULL, y = NULL, tol = sd(c(z)) * 0.1^2, max.sift = 20, 
    boundary = "reflexive", boundperc = 0.3, max.imf = 5, sm = "none", 
    smlevels = 1, spar = NULL, weight = NULL, plot.imf = FALSE) 
}
\arguments{
\item{z}{matrix of an image observed at (\code{x}, \code{y})}
\item{x, y}{locations of regular grid at which the values in \code{z} are measured}
\item{tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{boundary}{specifies boundary condition from ``none", ``symmetric" or ``reflexive".}
\item{boundperc}{expand an image by adding specified percentage of image at the boundary when boundary condition is 'symmetric' or 'reflexive'.}
\item{max.imf}{the maximum number of IMF's}
\item{sm}{specifies whether envelop is constructed by interpolation, thin-plate smoothing, Kriging, local polynomial smoothing, or loess. 
Use ``none" for interpolation, ``Tps" for thin-plate smoothing, ``mKrig" for Kriging, ``locfit" for local polynomial smoothing, or ``loess" for loess.
See Kim et al. (2012) for detalis.}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{spar}{specifies user-supplied smoothing parameter of thin-plate smoothing, Kriging, local polynomial smoothing, or loess.}
\item{weight}{deprecated.}
\item{plot.imf}{specifies whether each IMF is displayed. If \code{plot.imf=TRUE}, click the plotting area to start the next step.}
}
\details{
This function performs the bidimenasional empirical mode decomposition
utilizing extrema detection based
on the equivalence relation between neighboring pixels. See Kim et al. (2012) for detalis.
}
\value{ 
\item{imf}{two dimensional IMF's}
\item{residue}{residue image after extracting the IMF's}
\item{maxindex}{index of maxima}
\item{minindex}{index of minima}
\item{nimf}{number of IMF's} 
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Kim, D., Park, M. and Oh, H.-S. (2012) Bidimensional Statistical Empirical Mode Decomposition. \emph{IEEE Signal Prosessing Letters}, 
\bold{19}, 191--194, doi: 10.1109/LSP.2012.2186566.
}
\seealso{
\code{\link{extrema2dC}}, \code{\link{extractimf2d}}.
}
\examples{
data(lena)
z <- lena[seq(1, 512, by=4), seq(1, 512, by=4)]
image(z, main="Lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
#lenadecom <- emd2d(z, max.imf = 4)
#imageEMD(z=z, emdz=lenadecom, extrema=TRUE, col=gray(0:100/100))

### Test Image
ndata <- 128

x <- y <- seq(0, 9, length=ndata)
meanf1 <- outer(sin(2 * pi * x), sin(2 * pi * y))
meanf2 <- outer(sin(0.5 * pi * x), sin(0.5 * pi * y))
meanf <- meanf1 + meanf2

snr <- 2
set.seed(77)
zn <- meanf + matrix(rnorm(ndata^2, 0, sd(c(meanf))/snr), ncol=ndata)

rangezn <- range(c(meanf1, meanf2, meanf, zn))
par(mfrow=c(2,2), mar=0.1 + c(0, 0.25, 3, 0.25))
image(meanf1, main="high frequency component", xlab="", ylab="", zlim=rangezn, 
    col=gray(100:0/100), axes=FALSE)
image(meanf2, main="low frequency component", xlab="", ylab="", zlim=rangezn, 
    col=gray(100:0/100), axes=FALSE)
image(meanf, main="test image", xlab="", ylab="", zlim=rangezn, col=gray(100:0/100), axes=FALSE)
image(zn, main="noisy image", xlab="", ylab="", zlim=rangezn, col=gray(100:0/100), axes=FALSE)

#out <- emd2d(zn, max.imf=3, sm="locfit", smlevels=1, spar=0.004125)
#par(mfcol=c(3,1), mar=0.1 + c(0, 0.25, 0.25, 0.25)) 
#image(out$imf[[1]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)
#image(out$imf[[2]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)
#image(out$imf[[3]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)
}
\keyword{nonparametric}
