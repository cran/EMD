\name{hilbertspec}
\alias{hilbertspec}
\title{Hilbert Transform and Instantaneous Frequency}
\description{
This function calculates teh amplitude and instantaneous frequency
using Hilbert transforma.
}
\usage{
hilbertspec(xt, tt = NULL) 
}
\arguments{
\item{xt}{matrix of multiple signals. Each column represents a signal.}
\item{tt}{observation index or time index}
}
\details{
This function calculates teh amplitude and instantaneous frequency
using Hilbert transforma.
}
\value{
\item{amplitude}{matrix of amplitudes for multiple signals \code{xt}}
\item{instantfreq}{matrix of instantaneous frequencies for multiple signals \code{xt}}
\item{energy}{cumulative energy of multiple signals}
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.
}
\seealso{
\code{\link{spectrogram}}.
}
\examples{
tt <- seq(0, 0.1, length = 2001)[1:2000]           
f1 <- 1776; f2 <- 1000
xt <- sin(2*pi*f1*tt) * (tt <= 0.033 | tt >= 0.067) + sin(2*pi*f2*tt)
 
### Before treating intermittence
interm1 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE)  
### After treating intermittence
interm2 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE, 
interm=0.0007)

par(mfrow=c(2,1), mar=c(2,2,2,1))
test1 <- hilbertspec(interm1$imf)
spectrogram(test1$amplitude[,1], test1$instantfreq[,1])

test2 <- hilbertspec(interm2$imf, tt=tt)
spectrogram(test2$amplitude[,1], test2$instantfreq[,1])
}
\keyword{nonparametric}
