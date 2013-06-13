\name{qfa.epiplot}
\alias{qfa.epiplot}
\title{Makes an epistasis plot from the full results of qfa.epi}
\description{
Creates a scatterplot of control fitnesses on the x-axis and query fitnesses on the y-axis, with those deemed to be hits (by FDR adjusted p-value) coloured.  Essentially, this function assumes that the experiment consists of a series of paired fitness observations for a collection (typically genome-wide) of deletion mutations either single mutations (x-axis) or the same single mutation in combination with a common background or query mutation (y-axis).  Fitting a linear regression to all observations (forced through the origin) and searching for significant deviations from that regression is equivalent to searching for mutations which show significant deviation from a multiplicative model of genetic interaction.  Genes whose deletions deviate from this model significantly can be said to interact with the query gene.
}
\usage{
qfa.epiplot(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",
yylab="Query Fitness",mmain="Epistasis Plot",fmax=0)
}
\arguments{
  \item{results}{The results of interaction analysis returned by the qfa.epi function.}
  \item{qthresh}{The fdr adjusted cutoff point for determining hits.}
  \item{fitratio}{The ratio of background mutant fitness to wildtype fitness, from the genetic indepdence model.
If FALSE, this is estimated from results using linear regression.}
  \item{ref.orf}{ORF for a reference strain (typically wild-type or a surrogate), whose fitness will be marked on the control and query axes of the interaction plot (horizontal and vertical blue lines). HIS3 is the default strain.}
  \item{xxlab}{x axis label}
  \item{yylab}{y axis label}
  \item{mmain}{Plot label}
  \item{fmax}{Maximum fitness range for both x-axis (control axis) and y-axis (query axis).  If 0, axis ranges are automatically chosen to include all data points.}
}
\keyword{qfa}