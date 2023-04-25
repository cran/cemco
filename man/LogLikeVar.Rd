\name{LogLikeVar}
\alias{LogLikeVar}
\title{
Log likelihood of the CemCO algorithm with covariates effects on distributions means and distributions covariance matrices.
}
\description{
Returns the log-likelihood of the CemCO algorithm with covariates effects on distributions means and distributions covariance matrices.
}
\usage{
LogLikeVar(data, Y, phi, G, y_cov)
}
\arguments{
  \item{data}{
A numeric vector, matrix, or data frame of observations. Non-numerical values should be converted to integer or float (e.g. dummies). If matrix or data frame, rows and columns correspond to observations (n) and variables (P).
}
  \item{Y}{
numeric matrix of data to use as covariates. Non-numerical values should be converted to integer or float (e.g. dummies).
}
  \item{phi}{
list of fitted parameters in the same format as the output of the CemCO function.
}
  \item{G}{
An integer specifying the numbers of mixture components (clusters).
}
  \item{y_cov}{
numeric matrix of data to use as covariates for the covariance effect. Non-numerical values should be converted to integer or float (e.g. dummies).
}
}
\details{
Calculate the log likelihood of the mixture guassian distribution given by the CemCO algorithm with covariates effects on distributions means and distributions covariance matrices. This function is used in the optimization process of the EM algorithm used to estimate the CemCO parameters.
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
X = cbind(rnorm(10), rnorm(10))
Y = cbind(rnorm(10), rnorm(10))
K = 2

fit <- CemCOVar(X, Y, K, Y[,1], max_iter=2  , n_start=1, cores=1)
ll <- LogLikeVar(X, Y, fit[[1]], K, Y[,1])
}

