\name{semd}
\alias{semd}
\title{Statistical Empirical Mode Decomposition}
\description{
This function performs empirical mode decomposition using spline smoothing not interpolation for sifting process.
The smoothing parameter is automatically detemined by cross-validation.
}
\usage{
semd(xt, tt=NULL, cv.kfold, cv.tol=0.1^1, cv.maxiter=20, 
    emd.tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1", boundary="periodic", 
    smlevels=1, max.imf=10)
}
\arguments{
\item{xt}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{cv.kfold}{the number of fold of cross-validation}
\item{cv.tol}{tolerance for cross-validation}
\item{cv.maxiter}{maximum iteration for cross-validation}
\item{emd.tol}{tolerance for stopping rule of sifting. If \code{stoprule=type5}, the number of iteration for S stoppage criterion.}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting. The \code{type1} stopping rule indicates that absolute values of envelope mean must be less than the user-specified tolerance level 
in the sense that the local average of upper and lower envelope is zero.
The stopping rules \code{type2}, \code{type3}, \code{type4} and \code{type5} are the stopping rules given by equation (5.5) of Huang et al. (1998), equation (11a), 
equation (11b) and S stoppage of Huang and Wu (2008), respectively.} 
\item{boundary}{specifies boundary condition from ``none", ``wave", ``symmetric", ``periodic" or ``evenodd". See Zeng and He (2004) for \code{evenodd} boundary condition.}

\item{smlevels}{specifies which level of the IMF is obtained by smoothing spline.}
\item{max.imf}{the maximum number of IMF's}
}
\details{
This function performs empirical mode decomposition using spline smoothing not interpolation for sifting process.
The smoothing parameter is automatically detemined by cross-validation.
Optimization is done by golden section search. See Kim et al. (2012) for details.
}
\value{ 
\item{imf}{IMF's}
\item{residue}{residue signal after extracting IMF's from observations \code{xt}}
\item{nimf}{the number of IMF's} 
\item{optlambda}{smoothing parameter minimizing prediction errors of cross-validation} 
\item{lambdaconv}{a sequence of smoothing parameters for searching optimal smoothing papameter} 
\item{perr}{prediction errors of cross-validation according to \code{lambdaconv}} 
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
\code{\link{extractimf}}, \code{\link{emd}}.
}
\examples{
ndata <- 2048
tt <- seq(0, 9, length=ndata)                 
xt <- sin(pi * tt) + sin(2* pi * tt) + sin(6 * pi * tt)  + 0.5 * tt 
set.seed(1)
xt <- xt + rnorm(ndata, 0, sd(xt)/5)

\dontrun{
### Empirical Mode Decomposition by Interpolation
emdbyint <- emd(xt, tt, max.imf = 5, boundary = "wave")
### Empirical Mode Decomposition by Smoothing
emdbysm <- semd(xt, tt, cv.kfold=4, boundary="wave", smlevels=1, max.imf=5)

par(mfcol=c(6,2), mar=c(2,2,2,1), oma=c(0,0,2,0))                              
rangext <- range(xt); rangeimf <- rangext - mean(rangext)
plot(tt, xt, xlab="", ylab="", main="signal", ylim=rangext, type="l")
mtext("Decomposition by EMD", side = 3, line = 2, cex=0.85, font=2)
plot(tt, emdbyint$imf[,1], xlab="", ylab="", main="imf 1", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,2], xlab="", ylab="", main="imf 2", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,3], xlab="", ylab="", main="imf 3", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,4], xlab="", ylab="", main="imf 4", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,5]+emdbyint$residue, xlab="", ylab="", main="remaining signal",
    ylim=rangext, type="l")

plot(tt, xt, xlab="", ylab="", main="signal", ylim=rangext, type="l")
mtext("Decomposition by SEMD", side = 3, line = 2, cex=0.85, font=2)
plot(tt, emdbysm$imf[,1], xlab="", ylab="", main="noise", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,2], xlab="", ylab="", main="imf 1", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,3], xlab="", ylab="", main="imf 2", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,4], xlab="", ylab="", main="imf 3", ylim=rangeimf,  type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$residue, xlab="", ylab="", main="residue", ylim=rangext, type="l")}
}
\keyword{nonparametric}
