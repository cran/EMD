\name{extractimf2d}
\alias{extractimf2d}
\title{Bidimensional Intrinsic Mode Function}
\description{
This function extracts the bidimensional intrinsic mode function from given an image
utilizing extrema detection based
on the equivalence relation between neighboring pixels.
}
\usage{
extractimf2d(residue, x=NULL, y=NULL, nnrow=nrow(residue), 
    nncol=ncol(residue), tol=sd(c(residue))*0.1^2, 
    max.sift=20, boundary="reflexive", boundperc=0.3,  
    sm="none", spar=NULL, weight=NULL, check=FALSE)
}
\arguments{
\item{residue}{matrix of an image observed at (\code{x}, \code{y})}
\item{x, y}{locations of regular grid at which the values in \code{residue} are measured}
\item{nnrow}{the number of row of an input image}
\item{nncol}{the number of column of an input image}
\item{tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{boundary}{specifies boundary condition from ``none", ``symmetric" or ``reflexive".}
\item{boundperc}{expand an image by adding specified percentage of image at the boundary when boundary condition is 'symmetric' or 'reflexive'.}
\item{sm}{specifies whether envelop is constructed by interpolation, thin-plate smoothing, Kriging, local polynomial smoothing, or loess.
Use ``none" for interpolation, ``Tps" for thin-plate smoothing, ``mKrig" for Kriging, ``locfit" for local polynomial smoothing, or ``loess" for loess.}
\item{spar}{specifies user-supplied smoothing parameter of thin-plate smoothing, Kriging, local polynomial smoothing, or loess.}
\item{weight}{deprecated.}
\item{check}{specifies whether the sifting process is displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
}
\details{
This function extracts the bidimensional intrinsic mode function from given image
utilizing extrema detection based
on the equivalence relation between neighboring pixels. See Kim et al. (2012) for detalis.
See Kim et al. (2012) for detalis.
}
\value{ 
\item{imf}{two dimensional IMF}
\item{residue}{residue signal after extracting the finest IMF from \code{residue}}
\item{maxindex}{index of maxima}
\item{minindex}{index of minima}
\item{niter}{number of iteration obtaining the IMF} 
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Kim, D., Park, M. and Oh, H.-S. (2012) Bidimensional Statistical Empirical Mode Decomposition. \emph{IEEE Signal Processing Letters}, 
\bold{19}, 191--194, doi: 10.1109/LSP.2012.2186566.
}
\seealso{
\code{\link{extrema2dC}}, \code{\link{emd2d}}.
}
\examples{
data(lena)
z <- lena[seq(1, 512, by=4), seq(1, 512, by=4)]

\dontrun{
lenaimf1 <- extractimf2d(z, check=FALSE)}
}
\keyword{nonparametric}
