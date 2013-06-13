\name{logist}
\alias{logist}
\title{Logistic growth curve model}
\description{
Logistic Model(http://en.wikipedia.org/wiki/Logistic_function#In_ecology:_modeling_population_growth) for culture growth curves, describing variation in culture cell density with time t.
}
\usage{
logist(K,r,g,t)}
\arguments{
  \item{K}{Culture carrying capacity (AU).  Same units as (normalised) cell density observed in growth curve. }
  \item{r}{Culture growth rate parameter (per day).}
  \item{g}{Inoculum denisty (AU).  Same units as (normalised) cell density observed in growth curve.}
  \item{t}{Time since inoculation (d).}
}
\keyword{logist}