\name{emd.pred}
\alias{emd.pred}
\title{Prediction by EMD and VAR model}
\description{
This function calculates prediction values and confidence limits
using EMD and VAR (vector autoregressive) model.
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
using EMD and VAR (vector autoregressive) model. See Kim et al. (2008) for detalis.
}
\value{
\item{fcst}{prediction values}
\item{lower}{lower limits of prediction}
\item{upper}{upper limits of prediction}
}
\references{
Kim, D, Paek, S.-H. and Oh, H.-S. (2008) A Hilbert-Huang Transform Approach for Predicting Cyber-Attacks.
\emph{Journal of the Korean Statistical Society}, \bold{37}, 277--283, doi:10.1016/j.jkss.2008.02.006.
}
\keyword{nonparametric}
