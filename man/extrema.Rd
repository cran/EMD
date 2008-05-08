\name{extrema}
\alias{extrema}
\title{Finding Local Extrema and Zero-crossings}
\description{
This function indentifies extrema and zero-crossings.
}
\usage{
extrema(y, ndata = length(y), ndatam1 = ndata - 1)
}
\arguments{
\item{y}{input signal}
\item{ndata}{the number of observation}
\item{ndatam1}{the number of observation - 1}
}
\details{
This function indentifies extrema and zero-crossings.
}
\value{ 
\item{minindex}{matrix of time index at which local minima are attained. 
Each row specifies a starting and ending time index of a local minimum}
\item{maxindex}{matrix of time index at which local maxima are attained. 
Each row specifies a starting and ending time index of a local maximum.}
\item{nextreme}{the number of extrema} 
\item{cross}{matrix of time index of zero-crossings. 
Each row specifies a starting and ending time index of zero-crossings.} 
\item{ncross}{the number of zero-crossings} 
}
\seealso{
\code{\link{extrema2dC}}, \code{\link{extractimf}}, \code{\link{emd}}.
}
\examples{
y <- c(0, 1, 2, 1, -1, 1:4, 5, 6, 0, -4, -6, -5:5, -2:2)
#y <- c(0, 0, 0, 1, -1, 1:4, 4, 4, 0, 0, 0, -5:5, -2:2, 2, 2)
#y <- c(0, 0, 0, 1, -1, 1:4, 4, 4, 0, 0, 0, -5:5, -2:2, 0, 0)

plot(y, type = "b"); abline(h = 0)
extrema(y)
}
\keyword{nonparametric}
