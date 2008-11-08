\name{lena}
\docType{data}
\alias{lena}
\title{Gray Lena image}
\description{
A 512x512 gray image of Lena.
}
\usage{
data(lena)
}
\format{
A 512x512 matrix.
}
\examples{
data(lena)
image(lena, col=gray(0:100/100), axes=FALSE)
}
\keyword{datasets}
