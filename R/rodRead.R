################ Read and bind rod output into suitable form for likelihood #################
rod.read<-function(path=".",files=c(),inoctimes="BarcodeTimes.txt",background="",treatments=c(),barcodes=c(),
master.plates=c(),screen.names=c(),ORF2gene=""){
# Create environment with inoculation times; checks if path specified
if (length(strsplit(path,".")[[1]])>1){
	pathT<-1
	#inocenv<-bctimes(paste(path,inoctimes,sep="/"))
}else{
	#inocenv<-bctimes(inoctimes)
	pathT<-0
}
inocenv<-bctimes(inoctimes)
# If no files specified, use all .txt in working directory
if (length(files)==0){fs<-list.files(path=path,pattern=".txt")} else {fs<-files}
fs<-fs[(fs!='BarcodeTimes.txt')&(fs!='ORF2GENE.txt')&(fs!="model.txt")&(fs!="test.txt")]
print("List of data files to read:")
print(fs)
if (pathT==1){fs<-paste(path,fs,sep="/")}
bigd<-data.frame()
for (f in fs){
	print(paste("Reading",as.character(f)))
	# Find out which ROD version and print
	rodversion<-as.numeric(scan(f,nlines=1,what="character")[4])
	print(paste("ROD Export Version",rodversion))
	if (rodversion==5){
		# Read in tab deliminated rod output
		roddata<-read.delim(f,skip=1)
		# Only take specified treatments, or all if none specified
		# and do the same for other arguments
		if (length(treatments)>0){
			roddata<-roddata[(roddata$Treatments%in%treatments),]}
		# Only take specified barcodes, or all if none specified
		if (length(barcodes)>0){
			roddata<-roddata[(roddata$Barcode%in%barcodes),]}
		if (length(master.plates)>0){
			roddata<-roddata[(roddata$MasterPlate.Number%in%master.plates),]}
		if (length(screen.names)>0){
			roddata<-roddata[(roddata$Screen.Name%in%screen.names),]}
		# Add column for inoculation times
		roddata$Inoc.Time<-""; bcs<-unique(roddata$Barcode)
		for (b in bcs){
		roddata[roddata$Barcode==as.character(b),'Inoc.Time']<-get(as.character(b),envir=inocenv)}
		# Correct column names
		names(roddata)[3]<-"Row"; names(roddata)[5]<-"Col"; names(roddata)[14]<-"ORF"
		names(roddata)[4]<-"Growth"; names(roddata)[15]<-"Date.Time"
		# Bind 
		bigd<-rbind(bigd,roddata)} else {stop("Unkown ROD version")} #ifROD5
	} #fs
# Add column for genetic background
bigd$Background<-background
# Print checks so people are sure they're using the correct data #
print(paste("Number of Barcodes =",length(unique(bigd$Barcode))))
print("Treatments :")
print(unique(bigd$Treatments))
print("Media:")
print(unique(bigd$Medium))
print(paste("Number of ORFs =",length(unique(bigd$ORF))))
print("Screens:")
print(unique(bigd$Screen.Name))
platesize<-max(as.numeric(bigd$Row))*max(as.numeric(bigd$Col))
if (length(bigd$Date.Time)%%platesize!=0){
	warning("Number of photos not multiple of plate size")}
print("Inoculation DateTimes:")
print(unique(bigd$Inoc.Time))
# Add a column of times since inoculation
fmt="%Y-%m-%d_%H-%M-%S"
time<-as.POSIXlt(as.character(bigd$Date.Time),format=fmt)
inoc<-as.POSIXlt(as.character(bigd$Inoc.Time),format=fmt)
difft<-as.numeric(difftime(time,inoc,units="days"))
bigd$Expt.Time<-difft
# Add a column of gene names if conversion file available
if (ORF2gene!=""){
	gdict<-orf2gdict(ORF2gene)
	bigd$Gene<-sapply(as.character(bigd$ORF),orf2g,gdict)
}else{bigd$Gene<-bigd$ORF}
return(bigd)
}