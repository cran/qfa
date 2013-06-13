\name{qfa.plot}
\alias{qfa.plot}
\title{Plots fitted model and data for all the colonies in results of qfa.fit}
\description{
Produces a multipage pdf of growth curves.  Each page corresponds to a single plate and growth curves are arrayed on the page according to their position on the plate.  Both observations and the fitted growth curve from qfa.fit are shown for each culture.  Where available, various fitness estimates are displayed for each culture, together with culture genotype.  These .pdfs are useful for visually checking quality of model fit & data.
}
\usage{
qfa.plot(file,results,d,fmt="\%Y-\%m-\%d_\%H-\%M-\%S",barcodes=c(),
master.plates=c(),treatments=c(),screen.names=c(),screenIDs=c(),
maxg=0,maxt=0,logify=FALSE,densityCol="Growth",curves=TRUE,
ylabel="Cell density (AU)",ptype="p")}
\arguments{
  \item{file}{The file to output the plots to.}
  \item{results}{The output of qfa.fit which contains the fitted curve parameters of colony growth you wish to plot.}
  \item{d}{The original data.frame fed to qfa.fit containing all of the timecourse data}
  \item{fmt}{The format in which Date.Time of measurement and inoculation time are stored}
  \item{barcodes}{Plot only for the plates with barcodes in this character vector; all by default.}
  \item{master.plates}{Plot only for the plates from master.plates in this character vector; all by default.}
  \item{treatments}{Plot only for the plates with treatments in this character vector; all by default.}
  \item{screen.names}{Plot only for the plates with screen.names in this character vector; all by default.}
  \item{screenIDs}{Plot only for the plates with screenIDs in this character vector; all by default.}
  \item{maxg}{Upper cell density (y-axis limit) for all growth curve plots.  Default value is zero.  If maxg=0, then upper fitnesses are chosen automatically to show all datapoints.}
  \item{maxt}{Growth curve is plotted from time t = 0 to maxt.}
  \item{logify}{Boolean indicating whether growth curve plots should be on a semilog scale.  Cell density (y-axis) only.}
  \item{densityCol}{Name of column in data frame d which contains cell density estimate.  Note that the image analysis software Colonyzer provides several possible alternatives.}
  \item{curves}{Boolean indicating whether fitted model curves should be drawn.  Useful to set this to false when generating diagnostic, data-only growth curves.  To do this, set curves = FALSE and pass the same (data) object twice, as both the d and results arguments.  For example: qfa.plot("test.pdf",df,df,curves=FALSE)}
  \item{ylabel}{String for y-axis label on growth curve plots.}
  \item{ptype}{Plot type for data: "p": points, "l": lines, "b": both}
}
\keyword{qfa}