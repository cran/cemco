\name{LogLike}
\alias{LogLike}
\title{
Log likelihood of the CemCO algorithm with covariates effects on distributions means.
}
\description{
Returns the log-likelihood of the CemCO algorithm with covariates effects on distributions means.
}
\usage{
LogLike(data, Y, phi, G)
}
\arguments{
  \item{data}{
A numeric vector, matrix, or data frame of observations. Non-numerical values should be converted to integer or float (e.g. dummies). If matrix or data frame, rows and columns correspond to observations (n) and variables (P).
}
  \item{Y}{
numeric matrix of data to use as covariates. Non-numerical values should be converted to integer or float (e.g. dummies).
}
  \item{phi}{
list of fitted parameters in the same format as the output of the CemCO function
}
  \item{G}{
An integer specifying the numbers of mixture components (clusters)
}
}
\details{
Calculate the log likelihood of the mixture guassian with covariates effects on distributions means. This function is used in the optimization process of the EM algorithm used to estimate the CemCO parameters.
}
\value{
Return the value of the log likelihood.
}
\references{
Stage I non-small cell lung cancer stratification by using a model-based clustering algorithm with covariates, Relvas et al.
}
\author{
Relvas, C. & Fujita, A.
}
\examples{
set.seed(42)
X = cbind(rnorm(60), rnorm(60))
Y = cbind(rnorm(60), rnorm(60))
K = 2

fit <- CemCO(X, Y, K, max_iter=10, n_start=1, cores=1)
ll <- LogLike(X, Y, fit[[1]], K)

}
