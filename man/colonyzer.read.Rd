\name{colonyzer.read}
\alias{colonyzer.read}
\title{Reads raw cell density timecourse data from Colonyzer output files}
\description{
Reads in and binds together all of the Colonyzer output files in a directory, so they are checked and ready for bayesian or likelihood inference.
Colonyzer is an open source image analysis tool for quantifying cell densities on agar plates: http://research.ncl.ac.uk/colonyzer/
}
\usage{
colonyzer.read(path=".",files=c(),experiment="ExptDescription.txt",
ORF2gene="",libraries="LibraryDescriptions.csv",screenID="")
}
\arguments{
  \item{path}{The path to the folder containing the Colonyzer .dat files to be read: working directory by default.}
  \item{files}{Character vector giving locations of Colonyzer .dat files to be read - overrides path}
  \item{experiment}{Name of text file describing the inoculation times, library and plate number for unique plates.  Taken relative to path if specified.  File must be a tab-delimited text file with no header containing the following columns:
  \itemize{
\item Barcode - Plate identifier
\item Area - Culture area (pixels)
\item Spot Row - Plate row number for culture
\item Trimmed Area - Integrated Optical Density, sum of pixel intensities in culture area
\item Spot Column - Plate column number for culture
\item Intensity - Total pixel intensity for square tile containing culture
\item Edge Pixels - Number of pixels classified as culture on edge of square tile
\item Threshold - Pixel intensity threshold used for image segmentation (after lighting correction)
\item X Offset - x-oordinate of top left hand corner of culture tile (pixels)
\item Y Offset - y-coordinate of top left hand corner of culture tile (pixels)
\item Treatments - Temperature(s) at which cultures were grown
\item Medium - Nutrients in plate agar
\item Image Name - Full name at image capture (includes barcode and date-time)
\item ORF Name - Array deletion y-number
\item Date of Image	- Date of image capture
\item Screen Name - Name of screen (identifies biological repeats, and experiment)
\item Library Name	- Identifier for library configuration (identifies particular culture location)
\item MasterPlate Number - Library plate identifier
\item Timeseries order - Photograph number
\item Colony Color R - Culture red channel intensity
\item Colony Color G - Culture green channel intensity
\item Colony Color B - Culture blue channel intensity
\item Background Color R	- Background red channel intensity (for current tile)
\item Background Color G	- Background green channel intensity (for current tile)
\item Background Color B	- Background blue channel intensity (for current tile)
\item Edge length - Number of culture pixels classified as being microcolony edge pixels (useful for classifying contaminants)
\item Tile Dimensions X - Culture tile width (pixels)
\item Tile Dimensions Y - Culture tile height (pixels)
}
}
  \item{ORF2gene}{Path to a tab-delimited text file containing two columns (with no headers) associating unique, systematic strain identifiers (e.g. yeast ORF Y-numbers) with human readable gene names (e.g. standard names from SGD).}
  \item{libraries}{Tab-delimited text file describing the array strains present in each row-column coordinate of each plate in a series of rectangular arrayed libraries.  Header row format is: "Library	ORF 	Plate	Row	Column	Notes".  Columns are:
  \itemize{
  \item Library - Library identifier (e.g. SDLV1)
  \item ORF - Systematic strain identifier
  \item Plate - Plate number
  \item Row - Row number
  \item Column - Column number
  \item Notes - Optional strain notes (e.g. is strain especially sick or missing?)
  }
  }
  \item{screenID}{Unique experiment identifier (e.g. QFA00001)}
}
\value{
An R data.frame where each row corresponds to a single observation on a single colony, with the value of the growth measurement in 'Growth', and the date and time of the measurement in 'Date.Time'. Other information about the observation is stored in the other columns.  Several columns returned are direct copies of Colonyzer output and mapped as follows:
\itemize{
\item Image.Name - Image Name
\item Row - Spot Row
\item Col - Spot Column
\item X.Offset - X Offset
\item Y.Offset - Y Offset
\item Area - Area
\item Trimmed - Trimmed Area
\item Threshold - Threshold
\item Intensity - Intensity
\item Edge.Pixels - Edge Pixels
\item Colony.Color.R - Colony Color R
\item Colony.Color.G - Colony Color G
\item Colony.Color.B - Colony Color B
\item Background.Color.R - Background Color R
\item Background.Color.G - Background Color G
\item Background.Color.B - Background Color B
\item Edge.length - Edge length
\item Tile.Dimensions.X - Tile Dimensions X
\item Tile.Dimensions.Y - Tile Dimensions Y
}

Extra columns are automatically added as follows.  Some of this information is derived from auxiliary files passed to the function such as the experimental description file, the orf-gene dictionary and the library description file:
\itemize{
\item Growth - A cell density surrogate built from trimmed Area normalised by tile area and maximum achievable pixel intensity: Trimmed/(Tile.Dimensions.X*Tile.Dimensions.Y*255)
\item Barcode - Plate identifier, essentially image name with date time and file extension stripped
\item Date.Time - Date time of image capture in YYYY-MM-DD_hh-mm-ss format
\item Inoc.Time - Date time that plate was inoculated.  If plate is grown at a high temperature, date time at which plate was moved into high temperature incubator.  The assumption in this case being that negligible growth occurred before plate temperature was shifted the the target temperature.
\item Treatments - Treatments applied to plate (e.g. temperature)
\item Medium - Medium contained in agar (e.g. nutrients or drugs added to agar)
\item Screen.Name - Unique identifier for experiment (usually identifies repeat number also if multiple repeats carried out).
\item RepQuad - Identifier for experiments scaling down from 1536 format plates to 384, indicating which quadrant on the original 1536 source plate the current 384 format plate belongs to.
\item MasterPlate.Number - Identifies which plate in the source library (as described in the library description file) corresponds to the current plate
\item Timeseries.order - Ordinal describing which photograph captured
\item Library.Name - Identifies which of the libraries identified in the library description file was used to construct this plate
\item ORF - Unique systematic identifier for the genotype of the strain at this location (e.g. yeast Y-number), as defined by library description file
\item Gene - Standard, human readable genotype identifier for the strain at this location, as defined by the ORF-Gene dictionary
\item Background - Tag identifying experiment, typically used to construct file names and axes titles in plots
\item Expt.Time - Number of days passed between inoculation (start of experiment) and current time
}

Finally, as well as returning the object above, this function prints a small report to screen, summarising the data returned.  This includes number of unique barcodes read, number of photos read, number of genotypes in experiment, number of unique culture observations made, a list of treatments applied, a list of media used, a list of unique screen names (e.g. replicates carried out), the plate dimensions (e.g. 1536, 384 or 96 format) and a list of unique inoculation dates.
}
\keyword{qfa}