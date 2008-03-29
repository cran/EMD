\name{emd2d}
\alias{emd2d}
\title{Two dimensional Empirical Mode Decomposition}
\description{
This function performs two dimensional empirical mode decomposition.
}
\usage{
emd2d(z, x = NULL, y = NULL, tol = sd(c(z)) * 0.1^2, max.sift = 20, 
    boundary = "reflexive", boundperc = 0.3, max.imf = 5,  
    smlevels = 1, weight = 0, plot.imf = FALSE) 
}
\arguments{
\item{z}{matrix of an image observed at (\code{x}, \code{y})}
\item{x, y}{locations of regular grid at which the values in \code{z} are measured}
\item{tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{boundary}{specifies boundary condition}
\item{boundperc}{expand an image by adding specified percentage of image at the boundary when boundary condition is 'symmetric' or 'reflexive'.}
\item{max.imf}{the maximum number of IMF's}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{weight}{the smoothness of thin plate spline is determined by \code{weight} times smoothing parameter of GCV.}
\item{plot.imf}{specifies whether each IMF is displayed. If \code{plot.imf=TRUE}, click the plotting area to start the next step.}
}
\details{
This function performs two dimensional empirical mode decomposition.
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
}
\keyword{nonparametric}
