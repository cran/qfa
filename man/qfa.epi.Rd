\name{qfa.epi}
\alias{qfa.epi}
\title{Finds genetic interaction strengths and p-values}
\description{
Fits a genetic independence model between control strains and double mutant strains, either using rjags and a Bayesian linear regression model, or lm and maximum likelihood. For each ORF, the probability that it is a false discovery of a suppressor or enhancer is calculated. These probabilities are then fdr corrected and returned along with genetic interaction scores.
}
\usage{
qfa.epi(double,control,qthresh,orfdict="ORF2GENE.txt",
GISthresh=0.0,plot=TRUE,modcheck=TRUE,fitfunct=mdrmdp,wctest=TRUE)
}
\arguments{
  \item{double}{Either a qfa.posterior or the results of qfa.fit for the double mutants}
  \item{control}{Either a qfa.posterior or the results of qfa.fit for the control strains}
  \item{qthresh}{The FDR corrected cut off}
  \item{orfdict}{Location of file giving a column of ORFs first and a column of corresponding gene names second - so gene names can be plotted}
  \item{GISthresh}{When returning interaction hitlists, this variable determines the cutoff for strength of genetic interaction.}
  \item{plot}{If TRUE, then a 2-way fitness plot is made.}
  \item{modcheck}{If TRUE then diagnostic residual plots are output to \dQuote{ModelCheck.pdf}}
  \item{fitfunct}{The name of a fitness function whose arguments are, in order, (K,r,g,v) (carrying capacity, rate and initial size of population and shape parameter for generalised logistic growth model).}
  \item{wctest}{If TRUE, then use the Wilcoxon test for differences in medians as a measure of statistical significance of genetic interaction.  This is the default.  If FALSE, then use a t-test for difference in mean fitnesses instead.}
}
\value{
Returns an R list containing three data frames: Results, Enhancers and Suppressors.  Each data frame has the following columns:
\itemize{
	\item ORF - Unique strain genotype identifier (e.g. Y-number for yeast strains)
	\item Gene - Human readable genotype identifier
	\item P - p-value for significance of difference between control and query strain fitnesses
	\item Q - q-value for significance of difference between control and query strain fitnesses.  This is FDR corrected p-value
	\item GIS - Genetic interaction strength.  Deviation of (mean or median, depending on value of wctest) observed query strain fitness from expected fitness given control query strain fitness and a multiplicative model of genetic interaction.
	\item QueryFitnessSummary - Summary statistic for all available replicate observations of query strain fitness (mean or median, depending on value of wctest).
	\item ControlFitnessSummary - Summary statistic for all available replicate observations of control strain fitness (mean or median, depending on value of wctest).
	\item QuerySE - Standard error on mean of query strain fitness observations
	\item ControlSE - Standard error on mean of control strain fitness observations
	\item TestType - Type of statistical test for significant difference carried out (i.e. Wilcoxon or t-test)
	\item SummaryType - Type of summary statistic used for fitnesses (i.e. mean or median)
	\item cTreat - Treatment applied to control plates
	\item cMed - Medium added to agar in control plates
	\item cBack - Control plate background tag (experiment identifier)
	\item qTreat - Treatment applied to query plates
	\item qMed - Medium added to agar in query plates
	\item qBack - Query plate background tag (experiment identifier)
	\item Type - Type of genetic interaction observed (suppressor, enhancer, positive, negative).  This is assigned for strains with abs(GIS)>GISthresh and by comparing q-value with qthresh.
}
}
\keyword{qfa}