\name{pgis}
\alias{pgis}
\title{Calculate strength and significance of genetic interaction.}
\description{
This function carries out t-test (parametric) or Mann-Whitney test (non-parameteric) for the difference between a set of observed query strain fitnesses and the predicted query strain fitnesses, given a set of observed control strain fitnesses (see epistasis calculations in qfa.epi).  It deals with non-uniqueness of replicate observations (e.g. all replicate strains are dead, giving repeated values of zero for fitness) and insufficient numbers of replicates sensibly.  Too many tied observations can render tests invalid. 
}
\usage{
pgis(orf,m,cFs,dFs,wilcoxon)}
\arguments{
  \item{orf}{Open Reading Frame (ORF) of array library gene (typically a deletion) of interest.}
  \item{m}{Slope of linear regression through origin fit to all available (typically genome-wide) observations.}
  \item{cFs}{List of available replicate observations of control strain fitnesses, labelled by ORF.}
  \item{dFs}{List of available replicate observations of query (or double mutant) strain fitnesses, labelled by ORF.}
  \item{wilcoxon}{Boolean specifying whether to use the Mann-Whitney, non-parameteric test for significance of difference between control and query strains (TRUE) or the parametric t-test (FALSE).  Default is TRUE.}
}
\keyword{pgis}