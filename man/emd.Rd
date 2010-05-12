\name{emd}
\alias{emd}
\title{Empirical Mode Decomposition}
\description{
This function performs empirical mode decomposition.
}
\usage{
emd(xt, tt=NULL, tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1",  
    boundary="periodic", smlevels=c(1), sm="none", spar=NA, weight=20, 
    check=FALSE, max.imf=10, plot.imf=TRUE, interm=NULL)
}
\arguments{
\item{xt}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting}
\item{boundary}{specifies boundary condition from ``none", ``wave", ``symmetric", ``periodic" or ``evenodd".}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{sm}{specifies whether envelop is constructed by smoothing spline.} 
\item{spar}{specifies user-supplied smoothing parameter of spline.}
\item{weight}{the smoothness of spline is determined by \code{weight} times smoothing parameter of GCV.}
\item{check}{specifies whether the sifting process is displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
\item{max.imf}{the maximum number of IMF's}
\item{plot.imf}{specifies whether each IMF is displayed. If \code{plot.imf=TRUE}, click the plotting area to start the next step.}
\item{interm}{specifies vector of periods to be excluded from the IMF's.}
}
\details{
This function performs empirical mode decomposition.
}
\value{ 
\item{imf}{IMF's}
\item{residue}{residue signal after extracting IMF's from observations \code{xt}}
\item{nimf}{the number of IMF's} 
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.
}

\seealso{
\code{\link{extrema}}, \code{\link{extractimf}}.
}
\examples{
### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2

par(mfrow=c(3,1), mar=c(2,1,2,1))
try <- emd(xt2, tt2, boundary="wave")

### Ploting the IMF's
par(mfrow=c(3,1), mar=c(2,1,2,1))
X11(); par(mfrow=c(try$nimf+1, 1), mar=c(2,1,2,1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf) {
    plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
    main=paste(i, "-th IMF", sep="")); abline(h=0)
}
plot(tt2, try$residue, xlab="", ylab="", main="residue", type="l", 
axes=FALSE); box()
}
\keyword{nonparametric}
