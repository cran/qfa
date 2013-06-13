testSlope=function(x,y){
	if((length(x)>1)&(length(y)==length(x))){
		linfit=lm(y~1+x)
		return(as.numeric(linfit$coefficients[2]))
	}else{
		return(0)
	}
}

slopeFun=function(x,y,N){
	LDat=length(x)
	# Pad x and y with N NAs at either end to allow for edge effects
	x=c(rep(NA,N),x,rep(NA,N))
	y=c(rep(NA,N),y,rep(NA,N))
	smallList=function(i) testSlope(x[(i-N):(i+N)],y[(i-N):(i+N)])
	s=sapply((N+1):(LDat+N),smallList)
	return(s)
	
}

lin.fit<-function(tim,growth,inocguess=1,xybounds=1,inits=list(),logTransform=FALSE,dplot=FALSE,asking=FALSE,orf="",detectThresh=0,...){
#dplot=TRUE
#asking=TRUE

	x=tim
	y=growth
	x=x[y>0]
	y=y[y>0]
	y=log(y)
	s=slopeFun(x,y,1)
	Kest=max(growth)

	slopecutoff=0
	ldenscutoff=detectThresh
	#udenscutoff=0.07 
	udenscutoff=0.2*Kest
	
	xnew=x[(s>slopecutoff)&(y>log(ldenscutoff))&(y<log(udenscutoff))]
	ynew=y[(s>slopecutoff)&(y>log(ldenscutoff))&(y<log(udenscutoff))]
	snew=s[(s>slopecutoff)&(y>log(ldenscutoff))&(y<log(udenscutoff))]
	
	pars=c(Kest,0,inocguess,1,1)
	#try({
	if(length(xnew)>1){
		# Fit straight line to filtered data
		linfit=lm(ynew~1+xnew)
		lA0=as.numeric(linfit$coefficients[1])
		r=as.numeric(linfit$coefficients[2])
		A0=exp(lA0)
	}else{
		r=0; A0=inocguess; Kest=inocguess; lA0=log(A0)
	}
	if((is.na(A0))|(is.na(r))) {r=0; A0=inocguess; Kest=inocguess; lA0=log(A0)}
	pars=c(Kest,r,A0,1,1)

	if (dplot)	{
		par(ask = asking)
		dt=signif(24*log(2)/r,3)
		inoc=signif(A0,3)

		#orf=dat$ORF[1]
		op=par(mfrow=c(2,1))
		plot(x,y,xlab="Time (d)",ylab="log cell density (AU)",main="Density cutoffs (green)")
		points(xnew,ynew,col="red")
		#abline(h=log(ldenscutoff),col="green")
		abline(h=log(udenscutoff),col="green")
		plot(x,s,xlab="Time (d)",ylab="slope of log density (AU)",main="Slope cutoff (blue)")
		points(xnew,snew,col="red")
		abline(h=slopecutoff,col="blue")
		par(op)

		mlab=paste(orf,"doubling time:",dt,"h","inoculum:",inoc)
		plot(x,y,xlab="Time (d)",ylab="log cell density (AU)",main=mlab)
		abline(lA0,r,col="red")
		if(length(xnew)>1){
			plot(xnew,fitted(linfit),type="l",col="red",xlab="Time (d)",ylab="log cell density (AU)",main=mlab)
			points(xnew,ynew)
			plot(linfit,ask=asking)
		}
	}
	return(pars)
}
#interact=FALSE
##pdf("ExponentialFitLogScale.pdf")
#svg(filename = "ExponentialModel%03d.svg",
#    width = 7, height = 7, pointsize = 12,
#    onefile = FALSE, family = "sans", bg = "white",
#    antialias = c("default", "none", "gray", "subpixel"))
#
#query=read.delim("LIQGROWTH/ANALYSISOUT/LIQGROWTH_Raw.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
#dat=query[(query$Barcode=="K000092_030_001")&(query$Row==3)&(query$Col==4),]
#lin.fit(dat$Expt.Time,dat$Growth,dplot=TRUE,asking=interact)
#dat=query[(query$Barcode=="K000092_030_001")&(query$Row==3)&(query$Col==5),]
#lin.fit(dat$Expt.Time,dat$Growth,dplot=TRUE,asking=interact)
#dat=query[(query$Barcode=="K000092_030_001")&(query$Row==3)&(query$Col==6),]
#lin.fit(dat$Expt.Time,dat$Growth,dplot=TRUE,asking=interact)
#dat=query[(query$Barcode=="K000092_030_001")&(query$Row==3)&(query$Col==7),]
#lin.fit(dat$Expt.Time,dat$Growth,dplot=TRUE,asking=interact)
#dev.off()



