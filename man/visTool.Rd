\name{makeVisTool}
\alias{makeVisTool}
\title{Making the visualisation tool}
\description{
Function closure which generates a function which can be run to initiate a Dynamic, Interactive X-Y (DIXY) visualisation tool.
}
\usage{
makeVisTool()
}
\value{

Returns a function for generating an interactive plot to compare control and query fitnesses and visualise interaction strengths.  The function returned takes the following arguments:

  \item{groups}{Data frame describing the names and members of groups of genes to be highlighted together.  Contains three columns: "GroupName", "GroupID" and "GroupORFs" describing a label for the group of genes, an ID specifying the source of the group (e.g. manually added or taken from Benschop et al. (2010)) and a list of the systematic gene names (Y numbers) for members of the group.}
  \item{orf2gene}{Data frame decribing a one-to-one relationship between standard gene names and systematic gene names.  Contains two columns: "ORF" and "Gene".  This association is used to convert between standard and systematic gene names during searching.  Note that only one standard gene name is associated with each systematic name, and so, during any search for standard gene names, if you happen to choose a version which is not included in this list, your target gene will not be found.}
  \item{GISfiles}{A list of objects output from the qfa package function report.epi, each of which compares the fitnesses of a library of strains grown as part of a control experiment and as part of a query experiment.}

Interactive plots are generated for each of the report.epi files whose filenames are listed in the GISfiles argument.
   
DIXY instructions:

Windows mouse
~~~~~~~~~~~~~~~
L click: Highlight gene/Rotate text position
R click: SGD (or press 'w' on keyboard)
M click: Remove last gene (or press 'd' on keyboard)

Mac mouse
~~~~~~~~~~~~~~~
Click: Highlight gene/Rotate text position

Keyboard
~~~~~~~~~~~~~~~
Left/Right arrow: change plot
Up/Down arrow: change functional complex highlighted
u: add new genes to list of functional complexes
z: select tool (toggle on and off)
s: add selection 
c: clear selection
w: open last gene highlighted in SGD
d: unhighlight last gene highlighted
t: toggle colours indicating positive and negative interaction
r: begin zoom (now click on top left and bottom right of new zoomed plot)
p: print current plot to QFAVisualisation.pdf
q: quit
}
\examples{
\dontrun{
visTool=makeVisTool()
visToolDemo()}
}
\keyword{qfa}