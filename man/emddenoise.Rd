\name{emddenoise}
\alias{emddenoise}
\title{Denoising by EMD and Cross-Validation}
\description{
This function performs denoising by empirical mode decomposition and cross-validation.
}
\usage{
emddenoise(xt, tt = NULL, cv.index, cv.level, cv.tol = 0.1^3, 
    cv.maxiter = 20, by.imf = FALSE, emd.tol = sd(xt) * 0.1^2, 
    max.sift = 20, stoprule = "type1", boundary = "periodic", 
    smlevels = c(1), sm = "none", spar = NA, weight = 20, 
    check = FALSE, max.imf = 10, plot.imf = FALSE, interm = NULL)
}
\arguments{
\item{xt}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{cv.index}{test dataset index according to cross-validation scheme}
\item{cv.level}{levels to be thresholded}
\item{cv.tol}{tolerance for the optimization step of cross-validation}
\item{cv.maxiter}{maximum iteration for the optimization step of cross-validation}
\item{by.imf}{specifies whether shrinkage is performed by each IMS or not.} 
\item{emd.tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting}
\item{boundary}{specifies boundary condition}
\item{smlevels}{specifies which level of the IMF is obtained by smoothing other than interpolation.}
\item{sm}{specifies whether envelop is constructed by smoothing spline.} 
\item{spar}{specifies user-supplied smoothing parameter of spline.}
\item{weight}{the smoothness of spline is determined by \code{weight} times smoothing parameter of GCV.}
\item{check}{specifies whether the sifting process is displayed.}
\item{max.imf}{the maximum number of IMF's}
\item{plot.imf}{specifies whether each IMF is displayed.}
\item{interm}{specifies vector of periods to be excluded from the IMF's.}
}
\details{
This function performs denoising by empirical mode decomposition and cross-validation.
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
#try10 <- emddenoise(xt, cv.index=cv.index, cv.level=2, by.imf=TRUE)
#try10$optlambda
}
\keyword{nonparametric}
