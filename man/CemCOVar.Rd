\name{CemCOVar}
\alias{CemCOVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit CemCO algorithm with covariates effects on cluster centroids and covariance matrices.
}
\description{
Model-based clustering based on parameterized finite Gaussian mixture models with covariates effects on the distribution means and the distribution covariance matrices. Models are estimated by an EM algorithm.
}
\usage{
CemCOVar(data, y, G, y_cov, max_iter=100, n_start=20, cores=4)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{y_cov}{
numeric matrix of data to use as covariates for the covariance effect. Non-numerical values should be converted to integer or float (e.g. dummies)
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
This function optimizes the log likelihood of the CemCO algorithm with covariates effects on cluster centroids and covariance matrices using a implementation of the EM algorithm. The covariates associated with the distributions means and with the distributions covariance matrices do not need to be the same.
}
\value{
The function output is a list
\item{fitted parameters}{The estimated parameters of the CemCO algorithm, including clusters centroids, covariance matrix, covariate effects of each cluster and a priori probability of each cluster.}
\item{log likelihood}{The optimal log likelihood estimated by the model}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
X = cbind(rnorm(20), rnorm(20))
Y = cbind(rnorm(20), rnorm(20))
K = 2

\donttest{
fit <- CemCOVar(X, Y, K, Y[,1], max_iter=5, n_start=1, cores=1)
params <- fit[[1]]  ## fitted parameters
ll <- fit[[2]]  ## log likelihood
}
}

