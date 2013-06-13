showDemo<-function(demoname="telomereCap"){
	# Function to just display text in demo
	file.show(system.file("demo", paste(demoname,".R",sep=""), package = "qfa"))
}
	
fitnessReport<-function(treatment,outputfile,dataframe){
	# Summarises mean and median fitnesses for all orfs in an .fit object
	print(outputfile)
	fitdf=dataframe[dataframe$Treatment==treatment,]

	orflst=unique(as.character(fitdf$ORF))
	orflst=sort(orflst)

	#meanRes=by(fitdf[,27:31],fitdf$ORF,median)
	medRes=tapply(fitdf$fit,fitdf$ORF,median,na.rm=TRUE)
	meanRes=tapply(fitdf$fit,fitdf$ORF,mean,na.rm=TRUE)
	varRes=tapply(fitdf$fit,fitdf$ORF,var,na.rm=TRUE)
	numRes=tapply(fitdf$fit,fitdf$ORF,length)
	seRes=tapply(fitdf$fit,fitdf$ORF,sd,na.rm=TRUE)/sqrt(numRes)
	orflst=names(medRes)
	genelst=as.character(fitdf$Gene[match(orflst,as.character(fitdf$ORF))])

	results=data.frame(Gene=genelst,ORF=orflst,MedianFit=medRes,MeanFit=meanRes,VarianceFit=varRes,NumRepeats=numRes,SEFit=seRes)

	packs = data.frame(installed.packages(),stringsAsFactors=FALSE)
	vno=packs$Version[packs$Package=="qfa"]
	QFAversion=paste("R package version:",vno)
	rTreat=paste("Treatment:",paste(unique(fitdf$Treatment),collapse=" "))
	rMed=paste("Medium:",paste(unique(fitdf$Medium),collapse=" "))
	rScrID=paste("Screen ID:",paste(unique(fitdf$ScreenID),collapse=" "))
	rGen=paste("Screen name:",paste(unique(fitdf$Screen.Name),collapse=" "))
	rLib=paste("Libraries:",paste(unique(fitdf$Library.Name),collapse=" "))
	rClient=paste("Client:",paste(unique(fitdf$Client),collapse=" "))
	rUser=paste("User:",paste(unique(fitdf$User),collapse=" "))
	rDate=paste("Date:",paste(unique(fitdf$ExptDate),collapse=" "))
	spacer="#########################################################"
	header=c(QFAversion,
	rTreat,rMed,rScrID,rGen,rLib,rClient,rUser,rDate,
	spacer)
	write.table(header,outputfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(results,"tmp.txt",sep="\t",quote=FALSE,row.names=FALSE)
	file.append(outputfile,"tmp.txt")
	file.remove("tmp.txt")	
	return(results)
}

report.epi<-function(results,filename){
	# Eliminate spurious precision to make smaller files
	results[,3:9]=signif(results[,3:9],4)
	packs = data.frame(installed.packages(),stringsAsFactors=FALSE)
	vno=packs$Version[packs$Package=="qfa"]

	QFAversion=paste("R package version:",vno)
	sumType=paste("Summary type:",paste(unique(results$SummaryType),collapse=" "))
	testType=paste("Test type:",paste(unique(results$TestType),collapse=" "))
	cTreat=paste("Control treatment:",paste(unique(results$cTreat),collapse=" "))
	cMed=paste("Control medium:",paste(unique(results$cMed),collapse=" "))
	cScrID=paste("Control screen ID:",paste(unique(results$cScrID),collapse=" "))
	cGen=paste("Control screen name:",paste(unique(results$cGen),collapse=" "))
	cLib=paste("Control libraries:",paste(unique(results$cLib),collapse=" "))
	cClient=paste("Control client:",paste(unique(results$cClient),collapse=" "))
	cUser=paste("Control user:",paste(unique(results$cUser),collapse=" "))
	cDate=paste("Control date:",paste(unique(results$cDate),collapse=" "))
	qTreat=paste("Query treatment:",paste(unique(results$qTreat),collapse=" "))
	qMed=paste("Query medium:",paste(unique(results$qMed),collapse=" "))
	qScrID=paste("Query screen ID:",paste(unique(results$qScrID),collapse=" "))
	qGen=paste("Query screen name:",paste(unique(results$qGen),collapse=" "))
	qLib=paste("Query libraries:",paste(unique(results$qLib),collapse=" "))
	qClient=paste("Query client:",paste(unique(results$qClient),collapse=" "))
	qUser=paste("Query user:",paste(unique(results$qUser),collapse=" "))
	qDate=paste("Query date:",paste(unique(results$qDate),collapse=" "))
	spacer="#########################################################"
	header=c(QFAversion,sumType,testType,
	cTreat,cMed,cScrID,cGen,cLib,cClient,cUser,cDate,
	qTreat,qMed,qScrID,qGen,qLib,qClient,qUser,qDate,
	spacer)
	write.table(header,filename,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	results$SummaryType=results$testType=NULL
	results$cTreat=results$cMed=results$cScrID=results$cGen=results$cLib=results$cClient=results$cUser=results$cDate=NULL
	results$qTreat=results$qMed=results$qScrID=results$qGen=results$qLib=results$qClient=results$qUser=results$qDate=NULL
	write.table(results,"tmp.txt",sep="\t",quote=FALSE,row.names=FALSE)
	file.append(filename,"tmp.txt")
	file.remove("tmp.txt")	
}

correlationReport<-function(scrnms,dataframe,outputfile,aw=4,ah=4,fitmax=185){
        # Tests the correlation of all possible repeats of the same plate/treatment
        # Useful tool for searching for incorrect plate orientation, or misplaced/mislabelled plates
        # Can also give clues about plates with incorrect medium
        
        # This might fail for various reasons
        # But probably not a critical part of any workflow, so wrap in "try"
        try({
        # Generate correlation plot for every possible 2-way comparison
        poss=combn(scrnms,2)
        trts=unique(as.character(dataframe$Treatment))
        plates=unique(as.numeric(as.character(dataframe$MasterPlate.Number)))
	  plates=plates[order(plates)]
	  pnum=length(plates)
        # Adjust aw, ah to get all plates on one page (e.g. 4*4 gives 16 panels, good for 15 plate lib)
        pdummy=aw*ah-pnum
	  pdf(outputfile)
        corrs=c()
        op<-par(mfrow=c(aw,ah),cex.main=0.8,mar=c(1,1,1,1),mgp=c(0,0,0))

        for(t in trts){
        for(comb in 1:length(poss[1,])){
        for(p in plates){
                rep1=poss[1,comb]; rep2=poss[2,comb]
                r1=dataframe[(as.character(dataframe$Treatment)==t)&(as.character(dataframe$Screen.Name)==rep1)&(as.numeric(as.character(dataframe$MasterPlate.Number))==p),]
                r2=dataframe[(as.character(dataframe$Treatment)==t)&(as.character(dataframe$Screen.Name)==rep2)&(as.numeric(as.character(dataframe$MasterPlate.Number))==p),]
                r1=r1[order(r1$Screen.Name,r1$MasterPlate.Number,r1$Gene),]
                r2=r2[order(r2$Screen.Name,r2$MasterPlate.Number,r2$Gene),]
                
		    if((length(r1$fit)==length(r2$fit))&(length(r1$fit)>0)){
                cols=rainbow(max(as.numeric(dataframe$MasterPlate.Number),na.rm=TRUE))
                correlate=cor(r1$fit,r2$fit)
                corrs=rbind(corrs,c(rep1,rep2,t,p,correlate))
                ptitle=paste("Trt:",t,"Plate:",p,"Corr:",formatC(correlate,4))
                #print(paste(rep1,rep2,ptitle))
                plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=rep1,ylab=rep2,main=ptitle,axes=FALSE)
                abline(0,1,lwd=3,col="grey")
                #text(r1$fit,r2$fit,r1$Gene,col="black",pos=4,offset=0.1,cex=0.4)
                points(r1$fit,r2$fit,col=cols[r1$MasterPlate.Number],pch=16,cex=0.4)
                #print(c(sum(r1$Gene==r2$Gene),length(r1$Gene)))
		    }else{plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=rep1,ylab=rep2,main=paste("Trt:",t,"Plate:",p),axes=FALSE)}
        }
        if(pdummy>0) for(p in 1:pdummy) plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab="",ylab="",main="",axes=FALSE)
        }
        }
        par(op)

        corrs=as.data.frame(corrs,stringsAsFactors=FALSE)
        colnames(corrs)=c("rep1","rep2","trt","p","correlate")
        corrs$p=as.numeric(corrs$p)
        corrs$correlate=as.numeric(corrs$correlate)
        hist(corrs$correlate,xlim=c(0,1),xlab="Correlation Coefficient",ylab="Frequency",main=outputfile)
        })
	  dev.off()
}

plateBoxplots<-function(dataframe,outputfile,fitmax=185){
# Generates plate by plate boxplots of culture fitnesses
# Useful for identifying plates with medium problems
trts=unique(dataframe$Treatment)
scrns=unique(dataframe$Screen.Name)
trts=trts[order(trts)]
scrns=scrns[order(scrns)]
pdf(outputfile)
for (trt in trts){
	for (scr in scrns){
		dt=dataframe[(dataframe$Screen.Name==scr)&(dataframe$Treatment==trt),]
		dt=dt[order(as.numeric(as.character(dt$MasterPlate.Number))),]
		plateNums=unique(as.numeric(as.character(dt$MasterPlate.Number)))
		plateNums=as.character(plateNums)
		dt$MasterPlate.Number<-factor(dt$MasterPlate.Number,levels=plateNums)
		boxplot(dt$fit~dt$MasterPlate.Number,notch=TRUE,xlab="Plate Number",ylab="Fitness",col=rainbow(23),main=paste(trt,scr),ylim=c(0,fitmax),cex.axis=0.5)
}
}
dev.off()
}

