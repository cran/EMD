\name{extractimf}
\alias{extractimf}
\title{Intrinsic Mode Function}
\description{
This function extracts intrinsic mode function from given a signal.
}
\usage{
extractimf(residue, tt=NULL, tol=sd(residue)*0.1^2, max.sift=20, 
    stoprule="type1", boundary="periodic", sm="none", spar=NULL, 
    alpha=NULL, check=FALSE, weight=NULL)
}
\arguments{
\item{residue}{observation or signal observed at time \code{tt}}
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
\item{spar}{specifies user-supplied smoothing parameter of spline smoothing, kernel smoothing, or local polynomial smoothing.}
\item{alpha}{deprecated.}
\item{check}{specifies whether the sifting process is displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
\item{weight}{deprecated.}
}
\details{
This function extracts intrinsic mode function from given a signal.
}
\value{ 
\item{imf}{imf}
\item{residue}{residue signal after extracting the finest imf from \code{residue}}
\item{niter}{the number of iteration to obtain the \code{imf}} 
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
\code{\link{extrema}}, \code{\link{emd}}.
}
\examples{
### Generating a signal
ndata <- 3000
X11(); par(mfrow=c(1,1), mar=c(1,1,1,1))
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2
plot(tt2, xt2, xlab="", ylab="", type="l", axes=FALSE); box()

### Extracting the first IMF by sifting process
tryimf <- extractimf(xt2, tt2, check=FALSE)
}
\keyword{nonparametric}
