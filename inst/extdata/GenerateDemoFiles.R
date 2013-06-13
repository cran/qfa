sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }
sourceDir("F:\\QFA\\pkg\\qfa\\R")

# Assume working directory is LOGS
cfolder="URA3" # Control folder
qfolder="CDC13-1" # Query folder

# Open all cdc13-1 image analysis files

#### QUERY directories
rpath=paste(qfolder,"/IMAGELOGS",sep="")
rexp=paste(qfolder,"/AUXILIARY/ExptDescription.txt",sep="")
rORF=paste(qfolder,"/AUXILIARY/ORF2GENE.txt",sep="")
rlib=paste(qfolder,"/AUXILIARY/LibraryDescription.txt",sep="")
rback=qfolder
rout=paste(qfolder,"/ANALYSISOUT/",qfolder,"_Raw.txt",sep="")
rrod=paste(qfolder,"/RODOUTPUT/",qfolder,"_ROD.txt",sep="")

# Read in colonyzer files and associate with data from expt description etc.
query=colonyzer.read(path=rpath,experiment=rexp,ORF2gene=rORF,libraries=rlib,background=rback)
query=query[(query$MasterPlate.Number==15)&(query$Treatments==27),]
write.table(query[,1:19],"cdc13-1Query.dat",col.names=FALSE,quote=FALSE,sep="\t")

# Open all ura3D image analysis files

#### CONTROL directories
rpath=paste(cfolder,"/IMAGELOGS",sep="")
rexp=paste(cfolder,"/AUXILIARY/ExptDescription.txt",sep="")
rORF=paste(cfolder,"/AUXILIARY/ORF2GENE.txt",sep="")
rlib=paste(cfolder,"/AUXILIARY/LibraryDescription.txt",sep="")
rback=cfolder
rout=paste(cfolder,"/ANALYSISOUT/",cfolder,"_Raw.txt",sep="")
rrod=paste(cfolder,"/RODOUTPUT/",cfolder,"_ROD.txt",sep="")

# Read in colonyzer files and associate with data from expt description etc.
control=colonyzer.read(path=rpath,experiment=rexp,ORF2gene=rORF,libraries=rlib,background=rback)
control=control[(control$MasterPlate.Number==15)&(control$Treatments==27),]
write.table(control[,1:19],"URA3Control.dat",col.names=FALSE,quote=FALSE,sep="\t")


