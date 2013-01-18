\name{kospi200}
\docType{data}
\alias{kospi200}
\title{Korea Stock Price Index 200}
\description{
the weekly KOSPI 200 index from January, 1990 to February, 2007.
}
\usage{
data(kospi200)
}
\format{
A list of date and KOSPI200 index
}
\details{
See Kim and Oh (2009) for the analysis for kospi200 data using EMD.
}
\references{
Kim, D. and Oh, H.-S. (2009) A Multi-Resolution Approach to Non-Stationary Financial Time Series Using the Hilbert-Huang Transform.
\emph{The Korean Journal of Applied Statistics}, \bold{22}, 499--513. 
}

\examples{
data(kospi200)
names(kospi200)
plot(kospi200$date, kospi200$index, type="l")
}
\keyword{datasets}
