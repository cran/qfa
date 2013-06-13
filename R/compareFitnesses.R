# Comparing DT with MDP...
# Modified version of densCols
dens<-function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(blues9[-(1:3)])) 
{
    xy <- xy.coords(x, y)
    select <- is.finite(xy$x) & is.finite(xy$y)
    x <- cbind(xy$x, xy$y)[select, ]
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
    xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
    dens <- map$fhat[cbind(xbin, ybin)]
    dens[is.na(dens)] <- 0
    dens
}

mkPlt<-function(qf,cutoff=0.03,lab,...){
	# Discard "dead" cultures before calculating correlation, as it skews result
	qfc=signif(cor(qf$MDR[(qf$MDR>0)],qf$MDP[(qf$MDR>0)]),3)
	cramp=colorRampPalette(c("blue", "orange", "red"), space = "Lab")
	cols=densCols(qf$MDR,qf$MDP,nbin=256,bandwidth=0.1,colramp=cramp)
	densit=dens(qf$MDR,qf$MDP,nbin=256,bandwidth=0.1,colramp=cramp)
	# Discard text for high density points
	txt=qf$Gene
	txt[densit>cutoff]=""
	delt=rep(0,length(qf$Gene))
	getDelt<-function(MDR,MDP){
		medMDR=median(qf$MDR[(qf$MDP>0.90*MDP)&(qf$MDP<1.1*MDP)])
		if(MDR<medMDR){return(-1)}else{return(1)}		
	}
	delt[densit<=cutoff]=sapply(qf$MDR[densit<=cutoff],getDelt,qf$MDP[densit<=cutoff])
	
	# Move text to the left of cloud to left of points, same for right...
	#reg=lm(qf$MDP~qf$MDR)
	#A0=as.numeric(reg$coefficients[1])
	#m=as.numeric(reg$coefficients[2])
	qf$pos=2
	#qf$pos[qf$MDR>(qf$MDP-A0)/m]=4
	qf$pos[delt>0]=4
	plot(qf$MDR,qf$MDP,pch=16,cex=0.5,main=paste(lab,"Correlation:",qfc),col=cols,xlab="Doubling rate (MDR)",ylab="Doubling potential (MDP)",...)
	text(qf$MDR,qf$MDP,txt,cex=0.3,pos=qf$pos)
}

summarise<-function(data,wctest=TRUE){
	###### Get ORF median fitnesses for control & double #######
	print("Calculating median (or mean) fitness for each ORF")
	## LIK ##
	# Get orfs in question
	orfs<-as.character(data$ORF)
	orfs<-unique(orfs)
	getMDR<-function(orf,fitframe){
        orfd<-fitframe[fitframe$ORF==orf,]
        return(orfd$MDR)
	}
	getMDP<-function(orf,fitframe){
        orfd<-fitframe[fitframe$ORF==orf,]
        return(orfd$MDP)
	}
	getGene<-function(orf,fitframe){
        orfd<-fitframe[fitframe$ORF==orf,]
        return(orfd$Gene[1])
	}
	# Get lists with fitnesses for each repeat
	mdrs<-lapply(orfs,getMDR,data)
	names(mdrs)<-orfs
	mdps<-lapply(orfs,getMDP,data)
	names(mdps)<-orfs
	genes<-lapply(orfs,getGene,data)
	names(genes)<-orfs
	# Get means or medians for each ORF
	if(wctest){mdrsumm<-sapply(mdrs,median)}else{mdrsumm<-sapply(mdrs,mean)}
	names(mdrsumm)<-orfs
	if(wctest){mdpsumm<-sapply(mdps,median)}else{mdpsumm<-sapply(mdps,mean)}
	names(mdpsumm)<-orfs
	mdrSE<-sapply(mdrs,sterr)
	mdpSE<-sapply(mdps,sterr)
	res=data.frame(ORF=as.character(orfs),Gene=as.character(genes),MDR=as.numeric(mdrsumm),MDP=as.numeric(mdpsumm),MDRSE=as.numeric(mdrSE),MDPSE=as.numeric(mdpSE))
}
