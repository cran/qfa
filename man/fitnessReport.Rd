\name{fitnessReport}
\alias{fitnessReport}
\title{Fitness Report}
\description{
Summarises mean and median fitnesses for all orfs in qfa.fit object, and writes them to file with an appropriate header.
}
\arguments{
  \item{treatment}{Fitnesses are filtered by this value in the qfa.fit treatment column before summaries are generated.}
  \item{outputfile}{Output file name.}
  \item{dataframe}{Dataframe containing fitnesses to be summarised.  Typically output from qfa.fit function.  Must add a "fit" column before passing to this function.}
}
\value{
Writes summary fitnesses to file with header containing some metadata.
}
\keyword{qfa}