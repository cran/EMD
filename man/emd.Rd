\name{emd}
\alias{emd}
\title{Empirical Mode Decomposition}
\description{
This function performs empirical mode decomposition.
}
\usage{
emd(xt, tt=NULL, tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1",  
    boundary="periodic", sm="none", smlevels=c(1), spar=NULL, alpha=NULL, 
    check=FALSE, max.imf=10, plot.imf=FALSE, interm=NULL, weight=NULL)
}
\arguments{
\item{xt}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{tol}{tolerance for stopping rule of sifting. If \code{stoprule=type5}, the number of iteration for S stoppage criterion.}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting. The \code{type1} stopping rule indicates that absolute values of envelope mean must be less than the user-specified tolerance level 
in the sense that the local average of upper and lower envelope is zero.
The stopping rules \code{type2}, \code{type3}, \code{type4} and \code{type5} are the stopping rules given by equation (5.5) of Huang et al. (1998), equation (11a), 
equation (11b) and S stoppage of Huang and Wu (2008), respectively.} 
\item{boundary}{specifies boundary condition from ``none", ``wave", ``symmetric", ``periodic" or ``evenodd". See Zeng and He (2004) for \code{evenodd} boundary condition.}
\item{sm}{specifies whether envelop is constructed by interpolation, spline smoothing, kernel smoothing, or local polynomial smoothing. 
Use ``none" for interpolation, ``spline" for spline smoothing, ``kernel" for kernel smoothing, or ``locfit" for local polynomial smoothing.
See Kim et al. (2012) for detalis.
}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{spar}{specifies user-supplied smoothing parameter of spline smoothing, kernel smoothing, or local polynomial smoothing.}
\item{alpha}{deprecated.}
\item{check}{specifies whether the sifting process is displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
\item{max.imf}{the maximum number of IMF's}
\item{plot.imf}{specifies whether each IMF is displayed. If \code{plot.imf=TRUE}, click the plotting area to start the next step.}
\item{interm}{specifies vector of periods to be excluded from the IMF's to cope with mode mixing.}
\item{weight}{deprecated.}
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

Huang, N. E. and Wu, Z. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.

Kim, D., Kim, K.-O. and Oh, H.-S. (2012) Extending the Scope of Empirical Mode Decomposition using Smoothing.
\emph{EURASIP Journal on Advances in Signal Processing}, \bold{2012:168}, doi: 10.1186/1687-6180-2012-168.

Zeng, K and He, M.-X. (2004) A simple boundary process technique for empirical mode decomposition.
\emph{Proceedings of 2004 IEEE International Geoscience and Remote Sensing Symposium}, \bold{6}, 4258--4261.
}
\seealso{
\code{\link{extrema}}, \code{\link{extractimf}}.
}
\examples{
### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2

try <- emd(xt2, tt2, boundary="wave")

### Ploting the IMF's
par(mfrow=c(try$nimf+1, 1), mar=c(2,1,2,1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf) {
    plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
    main=paste(i, "-th IMF", sep="")); abline(h=0)
}
plot(tt2, try$residue, xlab="", ylab="", main="residue", type="l", axes=FALSE); box()
}
\keyword{nonparametric}
