\name{emd.pred}
\alias{emd.pred}
\title{Prediction by EMD and VAR model}
\description{
This function calculates prediction values and confidence limits
using EMD and VAR model.
}
\usage{
emd.pred(varpred, trendpred, ci = 0.95, figure = TRUE) 
}
\arguments{
\item{varpred}{prediction result of IMF's by VAR model.}
\item{trendpred}{prediction result of residue by polynomial regression model.}
\item{ci}{confidence interval level.}
\item{figure}{specifies whether prediction result is displayed.}
}
\details{
This function calculates prediction values and confidence limits
using EMD and VAR model.
}
\value{
\item{fcst}{prediction values}
\item{lower}{lower limits of prediction}
\item{upper}{upper limits of prediction}
}
\keyword{nonparametric}
