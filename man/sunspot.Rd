\name{sunspot}
\docType{data}
\alias{sunspot}
\title{Sunspot Data}
\description{
sunspot from 1610 through 1995. 
}
\usage{
data(sunspot)
}
\format{
A list of year and sunspot 
}

\references{
Oh, H-S, Ammann, C. M., Naveau, P., Nychka, D. and Otto-Bliesner, B. L. 
(2003) Multi-resolution time series analysis applied to solar irradiance
and climate reconstructions. 
\emph{Journal of Atmospheric and Solar-Terrestrial Physics}, \bold{65}, 191--201.
}

\examples{
data(sunspot)
names(sunspot)
plot(sunspot$year, sunspot$sunspot, type="l")
}
\keyword{datasets}
