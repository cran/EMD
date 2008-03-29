\name{cvtype}
\alias{cvtype}
\title{Generating test dataset index for cross-validation}
\description{
This function generates test dataset index for cross-validation.
}
\usage{
cvtype(n, cv.bsize=1, cv.kfold, cv.random=TRUE)
}
\arguments{
\item{n}{the number of observation}
\item{cv.bsize}{block size of cross-validation}
\item{cv.kfold}{the number of fold of cross-validation}
\item{cv.random}{whether or not random cross-validation scheme should be used. 
  Set \code{cv.random=TRUE} for random cross-validation scheme}
}
\details{
This function provides index of test dataset according to various cross-validation scheme.
One may construct K test datasets in a way that each testset consists of blocks of b
consecutive data. Set \code{cv.bsize = b} for this. 
To select each fold at random, set \code{cv.random = TRUE}.
}
\value{
matrix of which row is test dataset index for cross-validation.
}

\examples{
# Traditional 4-fold cross-validation for 100 observations
cvtype(n=100, cv.bsize=1, cv.kfold=4, cv.random=FALSE)
# Random 4-fold cross-validation with block size 2 for 100 observations
cvtype(n=100, cv.bsize=2, cv.kfold=4, cv.random=TRUE)
}
\keyword{nonparametric}
