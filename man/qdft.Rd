% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qfa4.1.R
\name{qdft}
\alias{qdft}
\title{Quantile Discrete Fourier Transform (QDFT)}
\usage{
qdft(y, tau, n.cores = 1, cl = NULL)
}
\arguments{
\item{y}{vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)}

\item{tau}{sequence of quantile levels in (0,1)}

\item{n.cores}{number of cores for parallel computing (default = 1)}

\item{cl}{pre-existing cluster for repeated parallel computing (default = \code{NULL})}
}
\value{
matrix or array of quantile discrete Fourier transform of \code{y}
}
\description{
This function computes quantile discrete Fourier transform (QDFT) for univariate or multivariate time series.
}
\examples{
y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
tau <- seq(0.1,0.9,0.05)
y.qdft <- qdft(y,tau)
# Make a cluster for repeated use
n.cores <- 2
cl <- parallel::makeCluster(n.cores)
parallel::clusterExport(cl, c("tqr.fit"))
doParallel::registerDoParallel(cl)
y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
y.qdft <- qdft(y1,tau,n.cores=n.cores,cl=cl)
y2 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
y.qdft <- qdft(y2,tau,n.cores=n.cores,cl=cl)
parallel::stopCluster(cl)
}
