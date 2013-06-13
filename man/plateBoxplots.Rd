\name{plateBoxplots}
\alias{plateBoxplots}
\title{Plate Boxplots}
\description{
Generates plate-by-plate boxplots of culture fitnesses.
Useful quality control step for identifying plates with potential medium problems.
}
\usage{
plateBoxplots(dataframe,outputfile,fitmax=185)
}
\arguments{
  \item{dataframe}{Dataframe containing fitnesses to be summarised.  Typically output from qfa.fit function.  Must add a "fit" column before passing to this function.}
  \item{outputfile}{Output file name.}
  \item{fitmax}{Upper limit for y-axis in plots (maximum observable fitness in this combination of experiment and fitness definition. }
}
\value{
Generates a multi-page .pdf report demonstrating within-plate culture variance for all plates, media, treatments and replicates in a QFA experiment.}
\keyword{qfa}