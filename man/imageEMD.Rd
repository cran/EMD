\name{imageEMD}
\alias{imageEMD}
\title{Plot of Two Dimensional Decomposition Result}
\description{
This function performs plotting of input image, IMF's, residue and extrema.
}
\usage{
imageEMD(z = z, emdz, extrema = FALSE, \dots)
}
\arguments{
\item{z}{matrix of an image}
\item{emdz}{decomposition result}
\item{extrema}{specifies whether the extrma is displayed according to the level of IMF}
\item{\dots}{the usual arguments to the image function}
}
\details{
This function performs plotting of input image, IMF's, residue and extrema.
}
\examples{
data(lena)
z <- lena[seq(1, 512, by=4), seq(1, 512, by=4)]
image(z, main="Lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
#lenadecom <- emd2d(z, max.imf = 4)
#imageEMD(z=z, emdz=lenadecom, extrema=TRUE, col=gray(0:100/100))
}
\keyword{nonparametric}
