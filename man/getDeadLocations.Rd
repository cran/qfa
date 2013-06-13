\name{getDeadLocations}
\alias{getDeadLocations}
\title{Find dead cultures in SGA plates (1536 format), and report their location in spotted plates (384 format).}
\description{
Given Colonyzer quantifications of SGA plates in 1536 format, locate dead cultures and report the 384-format locations of those cultures for stripping.}
\usage{
getDeadLocations(SGAFile,SGAExpt,CutoffFrac=0.0025)
}
\arguments{
  \item{SGAFile}{The Colonyzer output for the final SGA plates.  These are 1536 format plates, and the Colonyzer file should be a .dat tab delimited text file.}
  \item{SGAExpt}{Experiment description file for the SGA, linking barcode with plate number.}
  \item{CutoffFrac}{Optional argument for specifying the minimum  value for Growth (normalised IOD) corresponding to detection of cells on 1536 plate.}
}
\value{
A data frame containing columns ROW384, COLUMN384, REP384 (repeat or quadrant identifier) and PLATE which will be useful for stripping.
}
\keyword{qfa}