\name{cvimpute.by.mean}
\alias{cvimpute.by.mean}
\title{Imputation by the mean of the two adjacent values}
\description{
This function performs imputation by the mean of the two adjacent values for test dataset of cross-validation.
}
\usage{
cvimpute.by.mean(y, impute.index)
}
\arguments{
\item{y}{observation}
\item{impute.index}{test dataset index for cross-validation}
}
\details{
This function performs imputation by the mean of the two adjacent values for test dataset of cross-validation.
See Kim et al. (2012) for detalis.
}
\value{
\item{yimpute}{imputed values by the mean of the two adjacent values}
}
\references{
Kim, D., Kim, K.-O. and Oh, H.-S. (2012) Extending the Scope of Empirical Mode Decomposition using Smoothing.
\emph{EURASIP Journal on Advances in Signal Processing}, \bold{2012:168}, doi: 10.1186/1687-6180-2012-168.
}
\seealso{
\code{\link{cvtype}}, \code{\link{semd}}.
}
\keyword{nonparametric}
