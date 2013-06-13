# Huge function closure to keep visualisation tool environment out of global environment

makeVisTool=function(){
	globs<-new.env()

	# Plotting function
	epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",ymin=0,ymax=0,xmin=0,xmax=0,ecol="green",scol="red"){
		enhancers<-gethits(results,qthresh,type="E")
		suppressors<-gethits(results,qthresh,type="S")
		others<-results[(!results$ORF%in%enhancers$ORF)&(!results$ORF%in%suppressors$ORF),]
		plot(NULL,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
		xlab=xxlab,ylab=yylab,main=mmain,
		col=8,pch=19,cex=0.5)
		# Add line for genetic independence
		if (fitratio!=FALSE){abline(0,fitratio,lwd=2,col=8)} else {
			abline(0,lm.epi(results$QueryFitnessSummary,results$ControlFitnessSummary,modcheck=FALSE),lwd=2,col=8)}
		# Add 1:1 fitness line
		abline(0,1,lwd=2,lty=4,col=8)
		# Add reference ORF fitnesses lines
		if (ref.orf!=FALSE){
			reforf<-results[results$ORF==ref.orf,]
			abline(v=reforf$ControlFitnessSummary,col="lightblue",lwd=2)
			abline(h=reforf$QueryFitnessSummary,col="lightblue",lwd=2)}
		# Add points for non-suppressors & non-enhancers
		points(others$ControlFitnessSummary,others$QueryFitnessSummary,col="grey",cex=0.5,pch=19)
		#text(others$ControlFitnessSummary,others$QueryFitnessSummary,others$Gene,col="grey",pos=4,offset=0.1,cex=0.4)
		# Add suppressors & enhancers
		if(length(enhancers$ORF)>0){
			points(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,col=ecol,pch=19,cex=0.5)
			#text(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
		}
		if(length(suppressors$ORF)>0){
			points(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,col=scol,pch=19,cex=0.5)
			#text(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
		}
	}

	makePlot=function(datno,focusPlot=TRUE,...){
		if(globs$compno>0){
			compnm=paste(globs$compno,globs$GROUPS$GroupName[globs$compno],"\t",globs$GROUPS$GroupID[globs$compno])
		}else{
			compnm=""
		}
		maintitle=paste(globs$datlist[[datno]]$datname,compnm,sep="\n")
		# In case of overly-long medium description, remove medium from axis labels
		if((nchar(globs$datlist[[datno]]$cMed)<22)&(nchar(globs$datlist[[datno]]$qMed)<22)){
			xlab=paste(globs$datlist[[datno]]$cScrID,globs$datlist[[datno]]$cCli,globs$datlist[[datno]]$cScrNm,globs$datlist[[datno]]$cLibs,globs$datlist[[datno]]$cUse,globs$datlist[[datno]]$cDate,globs$datlist[[datno]]$cTreat,globs$datlist[[datno]]$cMed)
			ylab=paste(globs$datlist[[datno]]$qScrID,globs$datlist[[datno]]$qCli,globs$datlist[[datno]]$qScrNm,globs$datlist[[datno]]$qLibs,globs$datlist[[datno]]$qUse,globs$datlist[[datno]]$qDate,globs$datlist[[datno]]$qTreat,globs$datlist[[datno]]$qMed)
		}else{
			xlab=paste(globs$datlist[[datno]]$cScrID,globs$datlist[[datno]]$cCli,globs$datlist[[datno]]$cScrNm,globs$datlist[[datno]]$cLibs,globs$datlist[[datno]]$cUse,globs$datlist[[datno]]$cDate,globs$datlist[[datno]]$cTreat)
			ylab=paste(globs$datlist[[datno]]$qScrID,globs$datlist[[datno]]$qCli,globs$datlist[[datno]]$qScrNm,globs$datlist[[datno]]$qLibs,globs$datlist[[datno]]$qUse,globs$datlist[[datno]]$qDate,globs$datlist[[datno]]$qTreat)
		}
		globs$orftargs=globs$dat$ORF[globs$targs]
		globs$dat=globs$datlist[[datno]]$res
		globs$targs=match(globs$orftargs,globs$dat$ORF)
		epiplot(globs$dat,0.05,mmain=maintitle,xxlab=xlab,yylab=ylab,ymin=globs$fymin[datno],ymax=globs$fymax[datno],xmin=globs$fxmin[datno],xmax=globs$fxmax[datno],...)
		if(globs$compno>0){
			lst=strsplit(globs$GROUPS$GroupORFs[globs$compno]," ")[[1]]
			dcomp=globs$dat[globs$dat$ORF%in%lst,]
			if(length(dcomp$ORF)>0){
				points(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,col="purple",pch=16,cex=1)
				text(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,dcomp$Gene,pos=1,cex=0.75)
			}	
		}
		if(length(globs$targs)>0){
			points(globs$dat$ControlFitnessSummary[globs$targs],globs$dat$QueryFitnessSummary[globs$targs],col="blue",pch=16,cex=0.2)
			text(globs$dat$ControlFitnessSummary[globs$targs],globs$dat$QueryFitnessSummary[globs$targs],globs$dat$Gene[globs$targs],pos=globs$posits,cex=0.75)
		}
		if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
	}

	ratPlot=function(datno,focusPlot=TRUE,qthresh=0.05,ecol="green",scol="red"){
		if(globs$compno>0){
			compnm=paste(globs$compno,globs$GROUPS$GroupName[globs$compno],"\t",globs$GROUPS$GroupID[globs$compno])
		}else{
			compnm=""
		}
		maintitle=paste(globs$datlist[[datno]]$datname,compnm,sep="\n")
		xlab=paste(globs$datlist[[datno]]$cScrID)
		ylab=paste(globs$datlist[[datno]]$qScrID)
		axislab=paste("Fit ",ylab,"/Fit ",xlab,sep="") 
		
		globs$orftargs=globs$dat$ORF[globs$targs]
		if(globs$ind=="gindex") ratlab="Ranked by gene name (alphabetical order)"
		if(globs$ind=="oindex") ratlab="Ranked by Y-number (chromosome coordinate)"
		if(globs$ind=="cindex") ratlab="Ranked by control fitness"
		if(globs$ind=="qindex") ratlab="Ranked by query fitness"
		if(globs$ind=="rindex") ratlab="Ranked by query/control fitness ratio"
		globs$dat=globs$datlist[[datno]]$res
		globs$targs=match(globs$orftargs,globs$dat$ORF)
		plot(globs$dat[[globs$ind]],globs$dat$ratio,xlab=ratlab,log="y",ylab=axislab,main=maintitle,col="grey",cex=0.5,pch=19)
		if(globs$compno>0){
			lst=strsplit(globs$GROUPS$GroupORFs[globs$compno]," ")[[1]]
			globs$dcomp=globs$dat[globs$dat$ORF%in%lst,]
			if(length(globs$dcomp$ORF)>0){
				points(globs$dcomp[[globs$ind]],globs$dcomp$ratio,col="purple",pch=16,cex=1)
				text(globs$dcomp[[globs$ind]],globs$dcomp$ratio,globs$dcomp$Gene,pos=1,cex=0.75)
			}	
		}
		if(length(globs$targs)>0){
			points(globs$dat[[globs$ind]][globs$targs],globs$dat$ratio[globs$targs],col="blue",pch=16,cex=0.2)
			text(globs$dat[[globs$ind]][globs$targs],globs$dat$ratio[globs$targs],globs$dat$Gene[globs$targs],pos=globs$posits,cex=0.75)
		}
		if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
	}

	targFun=function(x,y,datobj){
		if(globs$ratioPlot){
			d=((x-datobj[[globs$ind]])/max(datobj[[globs$ind]]))^2+(y-datobj$ratio)^2
		}else{
			d=(x-datobj$ControlFitnessSummary)^2+(y-datobj$QueryFitnessSummary)^2
		}
		best=which.min(d)
		return(best)
	}

	mouse=function(buttons,x,y){
		if(0%in%buttons) {
			# Problem with Cairo coordinate system here...
			globs$xf=grconvertX(x,from="ndc",to="user")
			globs$yf=grconvertY(y,from="ndc",to="user")
			if(globs$zooming==TRUE){
				if(identical(globs$zoomTL,c(-99,-99))) {
					globs$zoomTL=c(globs$xf,globs$yf)
				}else{
					globs$zoomBR=c(globs$xf,globs$yf)
					# Zoom in!
					globs$fxmin[globs$datno]=globs$zoomTL[1]
					globs$fxmax[globs$datno]=globs$zoomBR[1]
					globs$fymin[globs$datno]=globs$zoomBR[2]
					globs$fymax[globs$datno]=globs$zoomTL[2]
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					# Reset zooming parameters
					globs$zooming=FALSE
					globs$zoomTL=c(-99,-99)
					globs$zoomBR=c(-99,-99)
				}
			}else{
				if(globs$sel==0){
					# Find nearest gene, highlight location, draw plot
					candidate=targFun(globs$xf,globs$yf,globs$dat)
					if(!candidate%in%globs$targs) {
						globs$targs=c(globs$targs,candidate)
						globs$posits=c(globs$posits,1)
					}else{
						globs$targ=which(globs$targs==candidate)
						globs$pnew=globs$posits[globs$targ]
						globs$posits[globs$targ]=(globs$pnew%%4)+1
					}
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					if(length(globs$selx)>0) drawSel(globs$selx,globs$sely)
				}else{
					globs$selx=c(globs$selx,globs$xf)
					globs$sely=c(globs$sely,globs$yf)
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					drawSel(globs$selx,globs$sely)
				}
			}
		}
		if(1%in%buttons) {
			# Delete the last gene which was highlighted
			globs$targs=globs$targs[-length(globs$targs)]
			globs$posits=globs$posits[-length(globs$posits)]
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(2%in%buttons) {
			# Open gene webpage on SGD
			globs$xf=grconvertX(x,from="ndc",to="user")
			globs$yf=grconvertY(y,from="ndc",to="user")
			candidate=targFun(globs$xf,globs$yf,globs$dat)
			browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",globs$dat$ORF[candidate],sep=""))
		}
	}

	keybd=function(key){
		if(key=="q") return(invisible(1))
		if(key=="c") {
			# Clear highlighted/selected points
			globs$targs=c()
			globs$posits=c()
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="d"){
			# Unhighlight the last gene which was highlighted
			globs$targs=globs$targs[-length(globs$targs)]
			globs$posits=globs$posits[-length(globs$posits)]
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="w"){
			globs$targ=globs$targs[length(globs$targs)]
			browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",globs$dat$ORF[globs$targ],sep=""))
		}
		if(key=="Right") {
			globs$datno<<-(globs$datno%%length(globs$datlist))+1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
			globs$orftargs=globs$dat$ORF[globs$targs]
			globs$dat=globs$datlist[[globs$datno]]$res
			globs$targs=match(globs$orftargs,globs$dat$ORF)
		}
		if(key=="Left") {
			globs$datno=((globs$datno-2)%%length(globs$datlist))+1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
			orftargs=globs$dat$ORF[globs$targs]
			globs$dat=globs$datlist[[globs$datno]]$res
			globs$targs=match(globs$orftargs,globs$dat$ORF)
		}
		if(key=="Up") {
			globs$compno=globs$compno+1
			if(globs$compno>length(globs$GROUPS$GroupORFs)) globs$compno=1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="Down") {
			globs$compno=globs$compno-1
			if(globs$compno<=0) globs$compno=length(globs$GROUPS$GroupORFs)
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="z") {
			if(globs$sel==0){
				globs$selx=c()
				globs$sely=c()
			}
			globs$sel=(globs$sel+1)%%3
			if(globs$sel==2){
				# Close polygon
				globs$selx=c(globs$selx,globs$selx[1])
				globs$sely=c(globs$sely,globs$sely[1])
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
				}
				drawSel(globs$selx,globs$sely)
				globs$sel=0
			}
		}
		if((key=="s")&(length(globs$selx)>0)){
			if(globs$ratioPlot){
				ptsind=point.in.polygon(globs$dat[[globs$ind]],globs$dat$ratio,globs$selx,globs$sely)
			}else{
				ptsind=point.in.polygon(globs$dat$ControlFitnessSummary,globs$dat$QueryFitnessSummary,globs$selx,globs$sely)
			}
			globs$tsel=which(ptsind==1)
			globs$targs=c(globs$targs,globs$tsel)
			globs$posits=c(globs$posits,rep(1,length(globs$tsel)))
			globs$selx=c()
			globs$sely=c()
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="p") {
			psname=sprintf("QFAVisualisation%04d.ps",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),psname),"\n")
			cairo_ps(psname)
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="n") {
			pngname=sprintf("QFAVisualisation%04d.png",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),pngname),"\n")
			png(pngname,width=480*2,height=480*2,pointsize=12*2)
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
			}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="m") {
			pdfname=sprintf("QFAVisualisation%04d.pdf",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),pdfname),"\n")
			pdf(pdfname)
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="t") {
			ENH=globs$Ecol
			SUP=globs$Scol
			globs$Scol=ENH
			globs$Ecol=SUP
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="r") {
		# Zoom mode (click on TL and BR)
			globs$zooming=TRUE
			globs$zoomTL=c(-99,-99)
			globs$zoomBR=c(-99,-99)
		}
		if(key=="u") {
		# Get user input for a new list of genes
			if (Sys.info()['sysname']=="Windows") bringToTop(which=-1)
			newgrp=getText(globs$ORFGENE)
			newdf=data.frame(GroupName=newgrp[["Label"]],GroupID="R VisTool")
			newdf$GroupORFs=paste(newgrp[["ORFs"]],collapse=" ")
			globs$GROUPS=rbind(newdf,globs$GROUPS)
			globs$compno=1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="l") {
		# Log ratio plot
		globs$ratioPlot=!globs$ratioPlot
		if(globs$ratioPlot){
			ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}else{
			makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}
		}
		if(key=="i") {
		# Change index for ratio plot
		globs$indNo=globs$indNo+1
		globs$ind=globs$indPoss[globs$indNo%%length(globs$indPoss)+1]
		if(globs$ratioPlot){ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)}
		}
		if(key=="b") {
		# Print data to console
		cat("*************\n")
		print(globs$datlist[[globs$datno]]$rept$V1)
		}
		return(NULL)
	}

	visTool<-function(groups,orf2gene,GISfiles){
		# Function for generating interactive plot
		globs$GROUPS=groups()
		globs$ORFGENE=orf2gene
		globs$datlist=list()
		
		for (f in 1:length(GISfiles)){	
			report=getResults(GISfiles[f])
			report$res$ratio=pmax(0.005,report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
			report$res$gindex=rank(as.character(report$res$Gene))
			report$res$oindex=rank(as.character(report$res$ORF))
			report$res$cindex=rank(report$res$ControlFitnessSummary)
			report$res$qindex=rank(report$res$QueryFitnessSummary)
			report$res$rindex=rank(report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
			report$datname="Mean Fitnesses (MDR * MDP), t-test"
			globs$datlist[[f]]=report
		}

		globs$fxmax=c();globs$fymax=c()
		globs$fxmin=c();globs$fymin=c()

		for(d in globs$datlist){
			globs$fxmax=c(globs$fxmax,d$fMax)
			globs$fymax=c(globs$fymax,d$fMax)
			globs$fxmin=c(0,globs$fxmin)
			globs$fymin=c(0,globs$fymin)
		}
		
		print(globs$fxmit)
			
		globs$n=0
		globs$targs=c()
		globs$posits=c()
		globs$sel=0
		globs$selx=c()
		globs$sely=c()
		globs$indPoss=c("gindex","oindex","rindex","cindex","qindex")
		globs$indNo=0
		globs$ind=globs$indPoss[globs$indNo%%length(globs$indPoss)+1]
		
		# Toggling between colour of suppressors and enhancers
		globs$Ecol="green"
		globs$Scol="red"

		cat("Windows mouse\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Left click: Highlight gene/Rotate text position\n")
		cat("Right click: Open SGD page for gene (alternatively, press 'w' on keyboard)\n")
		cat("Middle click: Remove last gene (alternatively, press 'd' on keyboard)\n\n")
		cat("Mac mouse\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Click: Highlight gene/Rotate text position\n\n")
		cat("Keyboard\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Left/Right arrow: Switch to next/previous fitness plot (comparison between different pair of experiments).\n")
		cat("Up/Down arrow: Change group of genes currently highlighted (in purple).\n")
		cat("u: Add new group of genes to list of highlightable groups.\n")
		cat("z: Toggle select tool on/off. Select genes for highlighting by drawing a polygon on plot.  Press 's' to add selection when finished.\n") 
		cat("s: Highlight genes encircled using select tool ('z').\n") 
		cat("c: Clear highlighting from currently selected genes.\n") 
		cat("w: Open SGD web-page for last gene highlighted.\n")
		cat("d: Remove highlighting from last gene highlighted.\n")
		cat("t: Toggle colours (red/green) indicating positive and negative interaction.\n")
		cat("l: Toggle plot style between fitness plot and log ratio plot.\n")
		cat("i: Toggle between log ratio plot styles.\n")
		cat("r: Enter zoom mode.  Click on top left and bottom right of area to inspect.\n")
		cat("p: Save current plot as vector graphic to QFAVisualisation.ps.  Other filetypes can also be generated - 'n': .png and 'm': .pdf\n")
		cat("b: Print experimental metadata report to console window\n")
		cat("q: Quit tool and print gene names currently selected to console window.\n")

		globs$datno=1
		globs$plotno=1
		globs$compno=0
		globs$pos=1
		globs$zooming=FALSE
		globs$zoomTL=c(-99,-99)
		globs$zoomBR=c(-99,-99)
		globs$ratioPlot=FALSE

		globs$dat=globs$datlist[[globs$datno]]$res

		# Check if running under windows, if not, force X11
		# Note that event handling doesn't work under Linux :(
		sysinf=Sys.info()
		#if( sysinf["sysname"]!="Windows") x11()
		#if( sysinf["sysname"]=="Linux") Cairo()
		x11()
		print(globs$datno)
		if(globs$ratioPlot){
			ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}else{
			makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}

		getGraphicsEvent(prompt="L click: Highlight/Rotate, R click: SGD, M click: Remove, Left/Right: Change plot, z: select tool, s: add selection, c: clear, q: quit", onMouseDown=mouse, onKeybd=keybd)
		cat("~~~~~~~~~~~~~~~\n")
		cat("\n")
		cat("List of genes currently selected:\n")
		cat(globs$dat$Gene[globs$targs])
		cat("\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("\n")
		cat("Systematic names for genes currently selected:\n")
		cat(globs$dat$ORF[globs$targs])
		cat("\n")
		cat("\n")
		cat("\n")
		#if( sysinf["sysname"]!="Windows") dev.off()
		dev.off()
	}
	
	return(visTool)
}

benschopFromSource<-function(){
	# Read in functionally related complexes
	# Largely from Benschopp Mol Cell 2010
	# Can add some complexes manually
	compfile=file.path(system.file(package = "qfa"),"extdata","FunctionalComplexes.txt",sep="")
	Benschopp=read.delim(compfile,stringsAsFactors=FALSE,sep="\t")
	Benschopp$CompList=gsub(";"," ",gsub("\\s","",Benschopp$Complex.members..systematic.name))
	res=data.frame(GroupName=Benschopp$X..Complex.name,GroupID="Func.")
	res$GroupORFs=Benschopp$CompList
	return(res)
}

buildBenschop<-function(){
	Benschop=read.delim(file.path(system.file(package = "qfa"),"extdata","Benschop.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(Benschop)
}

visToolDemo<-function(groupFun=buildBenschop){
	orfile=file.path(file.path(system.file(package = "qfa"),"extdata","ORF2GENE.txt.gz"))
	ORFGENE=read.delim(orfile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
	colnames(ORFGENE)=c("ORF","Gene")
	ORFGENE=ORFGENE[!duplicated(ORFGENE$ORF),]

	# Read in GIS files
	filenames=list.files(file.path(system.file(package = "qfa"),"extdata"),pattern="*GIS.txt.gz",full.names=TRUE)
	visTool=makeVisTool()
	visTool(groupFun,ORFGENE,filenames)
}

NAtoBlank=function(x){
	if(is.na(x)){
		return("")
	}else{
		return(x)
	}
}

getResults<-function(filename){
	res=read.delim(filename,skip=20,header=TRUE,stringsAsFactors=FALSE)
	hdr=read.delim(filename,nrows=19,header=FALSE,stringsAsFactors=FALSE)
	qfaVersion=NAtoBlank(strsplit(hdr[1,1],": ")[[1]][2])
	summType=NAtoBlank(strsplit(hdr[2,1],": ")[[1]][2])
	testType=NAtoBlank(strsplit(hdr[3,1],": ")[[1]][2])
	cTreat=NAtoBlank(strsplit(hdr[4,1],": ")[[1]][2])
	cMed=NAtoBlank(strsplit(hdr[5,1],": ")[[1]][2])
	cScrID=NAtoBlank(strsplit(hdr[6,1],": ")[[1]][2])
	cScrNm=NAtoBlank(strsplit(hdr[7,1],": ")[[1]][2])
	cLibs=NAtoBlank(strsplit(hdr[8,1],": ")[[1]][2])
	cCli=NAtoBlank(strsplit(hdr[9,1],": ")[[1]][2])
	cUse=NAtoBlank(strsplit(hdr[10,1],": ")[[1]][2])
	cDate=NAtoBlank(strsplit(hdr[11,1],": ")[[1]][2])
	qTreat=NAtoBlank(strsplit(hdr[12,1],": ")[[1]][2])
	qMed=NAtoBlank(strsplit(hdr[13,1],": ")[[1]][2])
	qScrID=NAtoBlank(strsplit(hdr[14,1],": ")[[1]][2])
	qScrNm=NAtoBlank(strsplit(hdr[15,1],": ")[[1]][2])
	qLibs=NAtoBlank(strsplit(hdr[16,1],": ")[[1]][2])
	qCli=NAtoBlank(strsplit(hdr[17,1],": ")[[1]][2])
	qUse=NAtoBlank(strsplit(hdr[18,1],": ")[[1]][2])
	qDate=NAtoBlank(strsplit(hdr[19,1],": ")[[1]][2])
	fMax=max(c(res$QueryFitnessSummary,res$ControlFitnessSummary))
	return(
	list(res=res,qfaVersion=qfaVersion,summType=summType,testType=testType,
	cTreat=cTreat,cMed=cMed,cScrID=cScrID,cScrNm=cScrNm,cLibs=cLibs,cCli=cCli,cUse=cUse,cDate=cDate,
	qTreat=qTreat,qMed=qMed,qScrID=qScrID,qScrNm=qScrNm,qLibs=qLibs,qCli=qCli,qUse=qUse,qDate=qDate,
	fMax=fMax,rept=hdr)
	)
}

buildGO<-function(){
	GO=read.delim(file.path(system.file(package = "qfa"),"extdata","GOAnnotation.txt.gz"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(GO)
}

getText=function(ORFGENE){
	x=list()
	cat("Input a list of gene names, separated by spaces & press enter:\n")
	usergenes=scan(what=character(),nlines=1)
	usergenes=toupper(usergenes)
	cat("Please input a label for your list of genes (no spaces) & press enter:\n")
	lab=scan(what=character(),nlines=1,n=1)
	orfs=c()
	genes=c()
	for(gene in usergenes){
		if(gene%in%ORFGENE$Gene){
			orfs=c(orfs,ORFGENE$ORF[ORFGENE$Gene==gene])
			genes=c(genes,gene)
		} else if (gene%in%ORFGENE$ORF) {
			orfs=c(orfs,gene)
			genes=c(genes,ORFGENE$Gene[ORFGENE$ORF==gene])			
		}
	}
	x[["Genes"]]=genes
	x[["ORFs"]]=orfs
	x[["Label"]]=lab
	print(x)
	return(x)
}

drawSel=function(selx,sely){
	if(length(selx)==1) {
		points(selx,sely,col="black",pch=16,cex=0.4)
	}else{
		points(selx,sely,col="black",lty=2,type="l")
	}
}	

