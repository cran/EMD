\name{extrema2dC}
\alias{extrema2dC}
\title{Finding Local Extrema}
\description{
This function finds the bidimensional local extrema based
on the equivalence relation between neighboring pixels. 
}
\usage{
extrema2dC(z, nnrow=nrow(z), nncol=ncol(z))
}
\arguments{
\item{z}{matrix of an input image}
\item{nnrow}{the number of row of an input image}
\item{nncol}{the number of column of an input image}
}
\details{
This function finds the bidimensional local extrema based
on the equivalence relation between neighboring pixels. See Kim et al. (2012) for detalis.
}
\value{ 
\item{minindex}{index of minima. Each row specifies index of local minimum.}
\item{maxindex}{index of maxima. Each row specifies index of local maximum.}
}
\references{
Kim, D., Park, M. and Oh, H.-S. (2012) Bidimensional Statistical Empirical Mode Decomposition. \emph{IEEE Signal Processing Letters}, 
\bold{19}, 191--194, doi: 10.1109/LSP.2012.2186566.
}
\seealso{
\code{\link{extrema}}, , \code{\link{extractimf2d}}, \code{\link{emd2d}}.
}
\examples{
data(lena)
z <- lena[seq(1, 512, by=4), seq(1, 512, by=4)]

par(mfrow=c(1,3), mar=c(0, 0.5, 2, 0.5))
image(z, main="Lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)    

example <- extrema2dC(z=z)
localmin <- matrix(256, 128, 128)
localmin[example$minindex] <- z[example$minindex]
image(localmin, main="Local minimum", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)

localmax <- matrix(0, 128, 128)
localmax[example$maxindex] <- z[example$maxindex]
image(localmax, main="Local maximum", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
}
\keyword{nonparametric}
