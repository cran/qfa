\name{normalisePlates}
\alias{normalisePlates}
\title{Normalising culture fitness by plate}
\description{
Sometimes estimated culture fitnesses vary systematically depending on the plate on which they are inoculated.  Agar in individual plates could come from different batches, and therefore have slightly different levels of nutrients or water.  Plates could be inoculated at different times, and stored at slightly different temperatures for example.  Depending on inoculation method, inoculation time specified may be less accurate for individual plates.  Any of these issues could effect simulated fitness slightly.  This function allows us to normalise culture fitnessses across plates to eliminate such effects.  It should only really be used for small differences.  The technical causes of larger differences should be corrected experimentally before analysis instead of attempting to normalise them away.

Starting with a data frame describing the output from the qfa.fit function (with optional added columns from the makeFitness function) makeFitness finds all unique treatments in that data frame, calculates a median value from the indicated column for all plates subjected to that treatment and then normalises the fitnesses of each culture on each plate so that the median fitness on each plate is equal to the median fitness for all plates which have undergone the given treatment.
}
\usage{
normalisePlates(d,column)}
\arguments{
  \item{d}{Dataframe (output from qfa.fit) for normalisation.}
  \item{column}{String name of column to normalise (typically the fitness measure of interest).  See qfa.fit and makeFitness help files for descriptions of available culture fitness measures.}
}
\keyword{normalisePlates}