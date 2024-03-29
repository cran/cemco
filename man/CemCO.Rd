\name{CemCO}
\alias{CemCO}

\title{

Fit CemCO algorithm using multiple threads of the machine
}
\description{
Model-based clustering based on parameterized finite Gaussian mixture models with covariates effects on the distribution means. Models are estimated by an EM algorithm running in multiple threads of the machine
}
\usage{
CemCO(data, y, G, max_iter=100, n_start=20, cores=4)
}

\arguments{
 \item{data}{
A numeric vector, matrix, or data frame of observations. Non-numerical values should be converted to integer or float (e.g. dummies). If matrix or data frame, rows and columns correspond to observations (n) and variables (P).
}
  \item{y}{
numeric matrix of data to use as covariates. Non-numerical values should be converted to integer or float (e.g. dummies).
}
  \item{G}{
An integer specifying the numbers of mixture components (clusters)
}
  \item{max_iter}{
maximum number of iterations of the EM optimization (default value equals to 100)
}
  \item{n_start}{
how many random sets should be chosen? (default value equals to 20)
}
  \item{cores}{
number of cores for EM optimization (default value equals to 4)
}
}
\details{
This function optimizes the log likelihood of the CemCO algorithm using a implementation of the EM algorithm. If categorial features need to be used, please create dummies or use another encode method.
}
\value{
The function output is a list
\item{fitted parameters}{The estimated parameters of the CemCO algorithm, including clusters centroids, covariance matrix, covariate effects of each cluster and a priori probability of each cluster.}
\item{log likelihood}{The optimal log likelihood estimated by the model}

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
params <- fit[[1]]  ## fitted parameters
ll <- fit[[2]]  ## log likelihood
}


