\name{emddenoise}
\alias{emddenoise}
\title{Denoising by EMD and Thresholding}
\description{
This function performs denoising by empirical mode decomposition and thresholding.
}
\usage{
emddenoise(xt, tt = NULL, cv.index, cv.level, cv.tol = 0.1^3, 
    cv.maxiter = 20, by.imf = FALSE, emd.tol = sd(xt) * 0.1^2, 
    max.sift = 20, stoprule = "type1", boundary = "periodic", 
    max.imf = 10)
}
\arguments{
\item{xt}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{cv.index}{test dataset index according to cross-validation scheme}
\item{cv.level}{levels to be thresholded}
\item{cv.tol}{tolerance for the optimization step of cross-validation}
\item{cv.maxiter}{maximum iteration for the optimization step of cross-validation}
\item{by.imf}{specifies whether shrinkage is performed by each IMS or not.} 
\item{emd.tol}{tolerance for stopping rule of sifting. If \code{stoprule=type5}, the number of iteration for S stoppage criterion.}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting. The \code{type1} stopping rule indicates that absolute values of envelope mean must be less than the user-specified tolerance level 
in the sense that the local average of upper and lower envelope is zero.
The stopping rules \code{type2}, \code{type3}, \code{type4} and \code{type5} are the stopping rules given by equation (5.5) of Huang et al. (1998), equation (11a), 
equation (11b) and S stoppage of Huang and Wu (2008), respectively.} 
\item{boundary}{specifies boundary condition from ``none", ``wave", ``symmetric", ``periodic" or ``evenodd". See Zeng and He (2004) for \code{evenodd} boundary condition.}
\item{max.imf}{the maximum number of IMF's}
}
\details{
This function performs denoising by empirical mode decomposition and cross-validation. See Kim and Oh (2006) for details.
}
\value{ 
\item{dxt}{denoised signal}
\item{optlambda}{threshold values by cross-validation}
\item{lambdaconv}{sequence of lambda's by cross-validation}
\item{perr}{sequence of prediction error by cross-validation}
\item{demd}{denoised IMF's and residue}
\item{niter}{the number of iteration for optimal threshold value}
}    
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Huang, N. E. and Wu, Z. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.

Kim, D. and Oh, H.-S. (2006) Hierarchical Smoothing Technique by Empirical Mode Decomposition (Korean). \emph{The Korean Journal of Applied Statistics}, \bold{19}, 319--330.

Zeng, K and He, M.-X. (2004) A simple boundary process technique for empirical mode decomposition.
\emph{Proceedings of 2004 IEEE International Geoscience and Remote Sensing Symposium}, \bold{6}, 4258--4261.
}
\seealso{
\code{\link{cvtype}}, \code{\link{emd}}.
}
\examples{
ndata <- 1024
tt <- seq(0, 9, length=ndata)
meanf <- (sin(pi*tt) + sin(2*pi*tt) + sin(6*pi*tt)) * (0.0<tt & tt<=3.0) + 
 (sin(pi*tt) + sin(6*pi*tt)) * (3.0<tt & tt<=6.0) +
 (sin(pi*tt) + sin(6*pi*tt) + sin(12*pi*tt)) * (6.0<tt & tt<=9.0)
snr <- 3.0
sigma <- c(sd(meanf[tt<=3]) / snr, sd(meanf[tt<=6 & tt>3]) / snr, 
sd(meanf[tt>6]) / snr)
set.seed(1)
error <- c(rnorm(sum(tt<=3), 0, sigma[1]), 
rnorm(sum(tt<=6 & tt>3), 0, sigma[2]), rnorm(sum(tt>6), 0, sigma[3]))
xt <- meanf + error 

cv.index <- cvtype(n=ndata, cv.kfold=2, cv.random=FALSE)$cv.index 

\dontrun{
try10 <- emddenoise(xt, cv.index=cv.index, cv.level=2, by.imf=TRUE)

par(mfrow=c(2, 1), mar=c(2, 1, 2, 1))
plot(xt, type="l", main="noisy signal")
lines(meanf, lty=2)
plot(try10$dxt, type="l", main="denoised signal")
lines(meanf, lty=2)}
}
\keyword{nonparametric}
