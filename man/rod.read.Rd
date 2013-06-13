\name{rod.read}
\alias{rod.read}
\title{Reading of ROD raw timecourse data.  Deprecated.}
\description{
Reads in and binds together all of the ROD output files in a directory, so they are checked and ready for bayesian or likelihood inference .  Deprecated.
}
\usage{
rod.read(path=".",files=c(),inoctimes="BarcodeTimes.txt",background="",
treatments=c(),barcodes=c(),master.plates=c(),screen.names=c(),ORF2gene = "")
}
\arguments{
  \item{path}{The path to the folder containing the ROD files to be read: working directory by default. Do not have other text files here.}
  \item{files}{Character vector giving locations of ROD files to be read - overrides path}
  \item{inoctimes}{A text file whose first column includes the barcodes in the ROD files and whose second column is the corresponding inoculation date.times. Taken relative to path if specified.}
  \item{background}{The genetic background of the colonies in the ROD files}
  \item{treatments}{Store data only for the plates with treatments in this character vecor; all by default.}
  \item{barcodes}{Store data only for the plates with barcodes in this character vecor; all by default.}
  \item{master.plates}{Store data only for the plates from master.plates in this character vecor; all by default.}
  \item{screen.names}{Store data only for the plates with screen.names in this character vecor; all by default.}
  \item{ORF2gene}{Path to a tab-delimited text file containing two columns (with no headers) associating unique, systematic strain identifiers (e.g. yeast ORF Y-numbers) with human readable gene names (e.g. standard names from SGD).}
}
\value{
An R data.frame where each row corresponds to a single observation on a single colony, with the value of the growth measurement in 'Growth', and the date and time of the measurement in 'Date.Time'. Other information about the observation is stored in the other columns.
}
\keyword{qfa}
