\name{report.epi}
\alias{report.epi}
\title{Normalising culture fitness by plate}
\description{
Outputs the results from the qfa.epi function to a tab-delimited text file together with a header describing the medium, treatment and background for the control and query experiment, whether mean or median fitness summaries are used, whether reported p-values for significance of genetic interaction strengths are the result of t-tests or Mann-Whitney tests and indicating which version of the R package was used to generate the results.
}
\usage{
report.epi(results,filename)}
\arguments{
  \item{results}{Dataframe describing genetic interactions (output from epi.fit) to be summarised in a text file.}
  \item{filename}{Path to file.}
}
\keyword{report.epi}