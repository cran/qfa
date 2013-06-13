\name{Glogist}
\alias{Glogist}
\title{Generalised Logistic growth curve model}
\description{
Calculates the Generalised Logistic Model (http://en.wikipedia.org/wiki/Generalised_logistic_function) version of culture doubling time as a function of time t.
}
\usage{
Glogist(K,r,g,v,t)}
\arguments{
  \item{K}{Culture carrying capacity (AU).  Same units as (normalised) cell density observed in growth curve. }
  \item{r}{Culture growth rate parameter (per day).}
  \item{g}{Inoculum denisty (AU).  Same units as (normalised) cell density observed in growth curve.}
  \item{v}{Shape parameter.  Recover logistic model with v = 1.}
  \item{t}{Time since inoculation (d).}
}
\keyword{Glogist}