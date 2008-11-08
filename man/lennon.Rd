\name{lennon}
\docType{data}
\alias{lennon}
\title{Gray John Lennon image}
\description{
A 256x256 gray image of John Lennon.
}
\usage{
data(lennon)
}
\format{
A 256x256 matrix.
}
\examples{
data(lennon)
image(lennon, col=gray(100:0/100), axes=FALSE)
}
\keyword{datasets}
