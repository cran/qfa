\name{correlationReport}
\alias{correlationReport}
\title{Correlation Report}
\description{
Generates all possible fitness correlation plots, on a plate-by-plate basis, comparing all possible combinations of replicates for a given medium, treatment and plate number.
Useful tool for searching for incorrect plate orientation, or misplaced/mislabelled plates.
Can also give clues about plates with incorrect medium.
}
\usage{
correlationReport(scrnms,dataframe,outputfile,aw=4,ah=4,fitmax=185)
}
\arguments{
  \item{dataframe}{Dataframe containing fitnesses to be summarised.  Typically output from qfa.fit function.  Must add a "fit" column before passing to this function.}
  \item{outputfile}{Output file name.}
  \item{scrnms}{Screen names to test for correlation problems}
  \item{fitmax}{Upper limit for both x and y axes in correlation plots.  Maximum observable fitness for this combination of experiments and fitness definitions.}
  \item{aw}{Number of horizontal panels per page in output report.}
  \item{ah}{Number of vertical panels per page in output report.}
}
\value{
Generates a multi-page .pdf report demonstrating two-way between-replicate correlation for all possible pairs of plates and for all media, treatments and replicates in a QFA experiment.}
\keyword{qfa}