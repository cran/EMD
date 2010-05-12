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
\examples{
data(kospi200)
names(kospi200)
plot(kospi200$date, kospi200$index, type="l")
}
\keyword{datasets}
