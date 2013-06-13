
getDeadLocations<-function(SGAFile,SGAExpt,CutoffFrac=0.0025){

	# Read in the SGA Colonyzer file
	sga=read.delim(SGAFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	colnames(sga)=c("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
	# Create barcode and date columns
	sga$BARCODE=substr(sga$FILENAME,1,11)
	sga$DATETIME=substr(sga$FILENAME,13,31)

      # Read in the SGA Experimental desciption file
	expt=read.delim(SGAExpt,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	# Using row.names speeds up indexing..?
	row.names(expt)=expt$Barcode
	getPlate<-function(barc) expt[barc,]$PlateNumber
	getLib<-function(barc) expt[barc,]$Library

	# Fill in the plate numbers and library names
	sga$PLATE=sapply(sga$BARCODE,getPlate)
	sga$LIBRARY=sapply(sga$BARCODE,getLib)
	# Fill in the repeat numbers for 384 format
	sga$REP384=((sga$COLUMN-1)%%2+1)+2*((sga$ROW-1)%%2)
	# Fill in the row and column numbers for 384 format
	sga$ROW384=ceiling(sga$ROW/2)
	sga$COLUMN384=ceiling(sga$COLUMN/2)

	# Now we can define what we mean by "dead" pinned SGA cultures
	# Get average tile dimensions
	xdim=mean(sga$XDIM)
	ydim=mean(sga$YDIM)
	minfit=CutoffFrac*xdim*ydim*255
	sgadead=sga[sga$TRIMMED<minfit,]
	return(sgadead)
}

#SGAFile="CDC13-1EXO1/IMAGELOGS/SGAFinal.dat"
#SGAExpt="CDC13-1EXO1/AUXILIARY/SGAExptDescription.txt"

#strip=getDeadLocations(SGAFile,SGAExpt)




