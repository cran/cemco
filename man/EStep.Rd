\name{EStep}
\alias{EStep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the E step of the CemCO algorithm with covariates effects on distributions means.
}
\description{
Implements the expectation step of EM algorithm for parameterized Gaussian mixture models with covariates effects on the distribution means. It is also used to calculate the a posteriori probability of each observation belong to each cluster.
}
\usage{
EStep(data, Y, phi, G)
}
%- maybe also 'usage' for other objects documented here.
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
Calculate the a posteriori probability of each observation belong to each cluster given the data and the current parameters estimation.
}
\value{
Returns a n x G numeric matrix where n represents the number of observations (number of rows of data) and G (the number of clusters). The value i, j represents the probability of the i-th observation belong to j-th cluster.
}
\references{
Stage I non-small cell lung cancer stratification by using a model-based clustering algorithm with covariates, Relvas et al.
}
\author{
Relvas, C. & Fujita, A.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(42)
X = cbind(rnorm(60), rnorm(60))
Y = cbind(rnorm(60), rnorm(60))
K = 2

fit <- CemCO(X, Y, K, max_iter=10, n_start=1, cores=1)
prob <- EStep(X, Y, fit[[1]], K)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

