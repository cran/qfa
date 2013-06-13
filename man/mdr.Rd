\name{mdr}
\alias{mdr}
\title{Maximum Doubling Rate (MDR) for Generalised Logistic Function}
\description{
Calculates the Generalised Logistic Model (http://en.wikipedia.org/wiki/Generalised_logistic_function) version of the Maximum Doubling Rate as presented in Addinall et al. 2011 (http://www.plosgenetics.org/doi/pgen.1001362).  MDR is the inverse of the doubling time at the beginning of the experiment, where competition between cultures is at its lowest level.
}
\usage{
mdr(K,r,g,v)}
\arguments{
  \item{K}{Culture carrying capacity (AU).  Same units as (normalised) cell density observed in growth curve. }
  \item{r}{Culture growth rate parameter (per day).}
  \item{g}{Inoculum denisty (AU).  Same units as (normalised) cell density observed in growth curve.}
  \item{v}{Shape parameter.  Recover logistic model with v = 0.}
}
\keyword{mdr}