\name{loapproxfun}
\alias{loapproxfun}
\title{Model free growth curve approximation}
\description{
This is a function closure.  Given a timeseries dataset (growth curve data) it returns an appropriate approximating function.  If a loess smoothing span parameter appropriate for the data capture frequency (frequency of photographs) is specified, the approximating function will be a smoothed version of the data in the range of observations.  For all points before the first observation, the approximating function takes the value of the first smoothed version of the data.  Simiarly, beyond the final observation, the function returns the smoothed version of the data at the final timepoint.  If an inappropriate span parameter is passed to this function it will return a linear interpolation approximating function instead.  This can be prefarable where the loess smoother would add spurious curves to datasets with sparse observations (e.g. data captured manually 2 or 3 times per day) and should give very similar results.
}
\usage{
loapproxfun(t,g,span)
}
\arguments{
  \item{t}{List of observation times.}
  \item{g}{List of cell density observations corresponding to the times in t.}
  \item{span}{Loess smoothing span.  If the user specifies too small a value for a given frequency of data capture, loess smoothing will not be possible and linear interpolation will be used instead.}
}
\value{
Returns a function of time t.
}
\examples{
t=c(0,1,2,3,4,5)
g=c(0,2,4,5,5,4)
# Span is too small, revert to linear interpolation
func1=loapproxfun(t,g,span=0.2)
curve(func1,0,5,xlab="Time",ylab="Cell density")
# Span is big enough
func2=loapproxfun(t,g,span=3)
curve(func2,0,5,col="red",add=TRUE)
points(t,g)
}
\keyword{loapproxfun}