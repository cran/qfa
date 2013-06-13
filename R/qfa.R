#### Define Phenotype for fitness ####
#mdr<-function(K,r,g,v) sapply((r*v)/(log((K^v-g^v)/(K^v-(2^v)*(g^v)))+log(2)*v),na2zero)
mdr<-function(K,r,g,v) sapply((r*v)/log(1-(2^v-1)/((2^v)*(g/K)^v-1)),na2zero)
mdp<-function(K,r,g,v) sapply(log(K/g)/log(2),na2zero)
mdrmdp<-function(K,r,g,v) mdr(K,r,g,v)*mdp(K,r,g,v)
# Doubling time for generalised logistic function, as a function of time t
dtl<-function(K,r,g,v,t) sapply((-(r*t*v) + log((2**v*(1 - (K/g)**v)*((1 + (-1 + (K/g)**v)/exp(r*t*v))**(-1/v))**v)/(-1 + 2**v*((1 + (-1 + (K/g)**v)/exp(r*t*v))**(-1/v))**v)))/(r*v),na2zero)

# Logistic growth model #
logist<-function(K,r,g,t) (K*g*exp(r*t))/(K+g*(exp(r*t)-1))

# Generalised logistic growth model #
Glogist<-function(K,r,g,v,t) K/(1+(-1+(K/g)^v)*exp(-r*v*t))^(1/v)

####### Normalise a column in a data frame
normalisePlates=function(d,column){
	d$Treatment=as.character(d$Treatment); d[,column]=as.numeric(d[,column]); d$Barcode=as.character(d$Barcode)
	for(trt in sort(unique(d$Treatment))){
		med=median(d[d$Treatment==trt,column])
		for (b in unique(d$Barcode[d$Treatment==trt])){
			datlst=as.numeric(d[(d$Barcode==b)&(d$Treatment==trt),column])
			p.med=median(datlst)
			if(p.med>=1) datlst=datlst*med/p.med # Avoid division by zero median fitness...
			d[(d$Barcode==b)&(d$Treatment==trt),column]=datlst
		}		
	}
	return(as.numeric(d[,column]))
}

####### Convert data datetime to time from start in days ########
tconv<-function(tstring,startt,fmt){
	t<-as.POSIXlt(as.character(tstring),format=fmt)
	return(as.numeric(difftime(t,startt,units="days")))
}

# Timing function (internal) #
trep<-function(thing,tobj){
	print(paste("Finished",thing,"in",round(tobj[3],2),"seconds"))
	print(paste(round(tobj[3]/60,2),"minutes"))
}

# Converts NAs to zeros # Should get rid of this - CONOR.
na2zero<-function(num) {if (is.na(num)) num<-0; return(num)}

###### Create inoculation time environment #####
bctimes<-function(bctimef){
	startframe=read.delim(bctimef,colClasses=c("character"))
	starts<-new.env(hash=TRUE)
	for (k in 1:length(startframe$Barcode)){
		st=startframe$Start.Time[k]
		assign(startframe$Barcode[k],st,envir=starts)
	}
	return(starts)
} # bctimes

# Convert barcode to start time #
bc2st<-function(bc,inocenv) get(as.character(bc),envir=inocenv)

### Create ORF2Gene dictionary ###
orf2gdict<-function(dictionary){
	cache=new.env(hash = TRUE)
	# Read in dictionary txt file 
	orf2g=read.delim(dictionary,header=FALSE,
	colClasses=c("character"),col.names=c("ORF","Gene"))
	orf2g$ORF=toupper(orf2g$ORF)
	z=apply(orf2g,1,orfun,cache)
	return(cache)
}

# Function called in orf2gdict #
orfun<-function(row,environ){
	orf<-as.character(row[1])
	gene<-as.character(row[2])
	return(assign(orf,gene,envir=environ))
} #orfun

### Function that converts orf 2 gene ###
orf2g<-function(orf,dictenv){get(orf,envir=dictenv)}

## Standard error of mean
sterr <- function(x) sqrt(var(x)/length(x))

################################################## Epistasis Function ###########################################################
qfa.epi<-function(double,control,qthresh,orfdict="ORF2GENE.txt",
	GISthresh=0.0,plot=TRUE,modcheck=TRUE,fitfunct=mdrmdp,wctest=TRUE){
	###### Get ORF median fitnesses for control & double #######
	print("Calculating median (or mean) fitness for each ORF")
	## LIK ##
	# Get orfs in question
	orfs<-as.character(double$ORF)
	orfs<-orfs[orfs%in%as.character(control$ORF)]
	orfs<-unique(orfs)
	# Get lists with fitnesses for each repeat
	#cFstats<-lapply(orfs,orfstat,control,fitfunct)
	cFstats<-lapply(orfs,orffit,control)
	names(cFstats)<-orfs
	#dFstats<-lapply(orfs,orfstat,double,fitfunct)
	dFstats<-lapply(orfs,orffit,double)
	names(dFstats)<-orfs
	# Get means or medians for each ORF
	if(wctest){cFms<-sapply(cFstats,median)}else{cFms<-sapply(cFstats,mean)}
	names(cFms)<-orfs
	if(wctest){dFms<-sapply(dFstats,median)}else{dFms<-sapply(dFstats,mean)}
	cSe<-sapply(dFstats,sterr)
	dSe<-sapply(cFstats,sterr)
	names(dFms)<-orfs
	### Fit genetic independence model ###
	m<-lm.epi(dFms,cFms,modcheck)
	print(paste("Ratio of background mutant fitness to wildtype fitness =",round(m,4)))
	###### Estimate probability of interaction #######
	print("Calculating interaction probabilities")
	pg<-sapply(orfs,pgis,m,cFstats,dFstats,wilcoxon=wctest)
	pg<-as.data.frame(t(pg))
	colnames(pg)=c("p","gis")
	p<-pg$p

	# Adjust for multiple comparisons
	q<-p.adjust(p,"fdr")
	# Get gene names if needed
	# lik
	genes<-as.character(double$Gene[match(orfs,double$ORF)]); 
	#genes<-genes[genes%in%as.character(control$Gene)]
	# If gene names missing, find them with orf2g
	if (is.null(genes)){# Create orf-gene dictionary
		orfdict<-orf2gdict(orfdict)
		genes<-sapply(orfs,orf2g,orfdict)
	}
	
	# Get genetic interaction scores
	#meandiff<-mean(dFms-cFms)
	meandiff=1
	#gis<-dFms/mean(dFms)-cFms/mean(cFms)
	gis<-pg$gis/meandiff
	# Put into data.frame
	nObs=length(p)
	if(wctest) {testType="wilcoxon"; sumType="median"}else{testType="t-test"; sumType="mean"}
	results<-data.frame(ORF=orfs,Gene=genes,P=p,Q=q,GIS=gis,
	QueryFitnessSummary=dFms,ControlFitnessSummary=cFms,QuerySE=dSe,ControlSE=cSe,
	TestType=rep(testType,nObs),SummaryType=rep(sumType,nObs),
	cTreat=rep(control$Treatment[1],nObs),cMed=rep(control$Medium[1],nObs),cScrID=rep(control$ScreenID[1],nObs),
	qTreat=rep(double$Treatment[1],nObs),qMed=rep(double$Medium[1],nObs),qScrID=rep(double$ScreenID[1],nObs),
	qGen=rep(double$Screen.Name[1],nObs),cGen=rep(control$Screen.Name[1],nObs),
	cLib=rep(paste(unique(control$Library.Name),collapse="_"),nObs),qLib=rep(paste(unique(double$Library.Name),collapse="_"),nObs))
	if("Client"%in%colnames(control)) results$cClient=rep(control$Client[1],nObs)
	if("Client"%in%colnames(double)) results$qClient=rep(double$Client[1],nObs)
	if("ExptDate"%in%colnames(control)) results$cDate=rep(control$ExptDate[1],nObs)
	if("ExptDate"%in%colnames(double)) results$qDate=rep(double$ExptDate[1],nObs)
	if("User"%in%colnames(control)) results$cUser=rep(control$User[1],nObs)
	if("User"%in%colnames(double)) results$qUser=rep(double$User[1],nObs)
	
	results$Type<-apply(results,1,typemake,m)
	results<-results[order(results$GIS,results$Q,results$Type),]
	# Get rid of duplicate entries in results
	orflist<-unique(as.character(results$ORF))
	results<-results[match(orflist,results$ORF),]
	# Plot results
	if (plot==TRUE){qfa.epiplot(results,qthresh,m)}
	final<-list(Results=results,
	Enhancers=gethits(results,qthresh,type="E",GISthresh=GISthresh),
	Suppressors=gethits(results,qthresh,type="S",GISthresh=GISthresh))
	return(final)
}

############### Epistasis Functions ##################
# Makes epistasis plot for a given fdr level #
qfa.epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",fmax=0){
	enhancers<-gethits(results,qthresh,type="E")
	suppressors<-gethits(results,qthresh,type="S")
	others<-results[(!results$ORF%in%enhancers$ORF)&(!results$ORF%in%suppressors$ORF),]
	if (fmax==0){
		# Get plot parameters
		ymax=1.1*max(results$QueryFitnessSummary); ymin=0
		xmax=1.1*max(results$ControlFitnessSummary); xmin=0
	}else{
		ymin=0;ymax=fmax
		xmin=0;xmax=fmax
	}
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
	text(others$ControlFitnessSummary,others$QueryFitnessSummary,others$Gene,col="grey",pos=4,offset=0.1,cex=0.4)
	# Add suppressors & enhancers
	if(length(enhancers$ORF)>0){
		points(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,col='green',pch=19,cex=0.5)
		text(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
	if(length(suppressors$ORF)>0){
		points(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,col='red',pch=19,cex=0.5)
		text(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
}

## Extract hits from epistasis results object ##
gethits<-function(results,qthresh,type="S",all.types=FALSE,GISthresh=0){
	results<-results[results$Q<qthresh,]
	if (all.types==FALSE){results<-results[results$Type==type,]}
	results<-results[abs(results$GIS)>=GISthresh,]
	results<-results[!is.na(results$ORF),]
	return(results[order(results$GIS,results$Q,results$Type),])
}

# Function to calculate fitnesses for repeats
orfstat<-function(orf,fitframe,fitfunct){
	orfd<-fitframe[fitframe$ORF==orf,]
	return(fitfunct(orfd$K,orfd$r,orfd$g,orfd$v))
}

# Function to extract fitnesses from data frame with "fit" column
orffit<-function(orf,fitframe){
	orfd<-fitframe[fitframe$ORF==orf,]
	return(orfd$fit)
}

# Linear regression genetic ind. model fit
lm.epi<-function(doubles,controls,modcheck){
	#indmod<-lm(doubles~controls)
	indmod<-lm(doubles~0+controls)
	m<-as.numeric(indmod$coefficients)
	# Check if linear regression OK
	if (modcheck==TRUE){
		pdf("ModelCheck.pdf")
		resids<-indmod$residuals
		hist(resids,xlab="Fitness Residuals",main="Histogram of Residuals")
		qqnorm(resids/sd(resids),main="Normal QQ Plot")
		abline(0,1,lwd=2,col="red")
		dev.off()}
	return(m)
}

# Estimates p-value and estimated strength of interaction
pgis<-function(orf,m,cFs,dFs,wilcoxon=TRUE){
	# If this orf is not present in both lists, return appropriate p,gis
	if((length(dFs[[orf]])==0)|(length(cFs[[orf]])==0)){return(c(1,0))}
	# If fitnesses are not unique in one condition (e.g. all dead) and only one repeat in another (e.g. after stripping)
	# This would cause wilcoxon test to fail due to ties, and t.test to fail due to insufficient y?
	ldFS=length(dFs[[orf]]); lcFS=length(cFs[[orf]])
	ludFS=length(unique(dFs[[orf]])); lucFS=length(unique(cFs[[orf]]))
	if (((ludFS==1)&(ldFS>1)&(lcFS==1))|((lucFS==1)&(lcFS>1)&(ldFS==1))){return(c(1,median(dFs[[orf]],na.rm=TRUE)-m*median(cFs[[orf]],na.rm=TRUE)))}
	if ((ludFS==1)&(lucFS==1)){return(c(1,median(dFs[[orf]],na.rm=TRUE)-m*median(cFs[[orf]],na.rm=TRUE)))}
	# If both sets of cultures are dead (and fitnesses therefore equal) return appropriate p,gis
	if(sum(dFs[[orf]]==m*cFs[[orf]])==length(dFs[[orf]]==m*cFs[[orf]])){
		return(c(1,0))
	}else{
		if (wilcoxon){
			# Returns p-value for significance of difference, and estimate of difference between medians
			ctest<-wilcox.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
			p<-as.numeric(ctest$p.value)
			diff<-as.numeric(ctest$estimate)
			return(c(p,diff))
		}else{
			# t-test fails if only one element in either list, do one sample test
			if((ldFS<=1)|(lcFS<=1)){
				ctest<-t.test(dFs[[orf]]-m*cFs[[orf]])
				p<-as.numeric(ctest$p.value)
				diff<-as.numeric(ctest$estimate)
				return(c(p,diff))
			}else{
				# Returns p-value for significance of difference, and the difference between the means
				ctest<-t.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
				p<-as.numeric(ctest$p.value)
				diff<-as.numeric(ctest$estimate)
				diff<-diff[1]-diff[2]
				return(c(p,diff))
			}
		}
	}
}

# Get type of interaction
typemake<-function(row,m){
	summd<-as.numeric(row['QueryFitnessSummary'])
	summc<-as.numeric(row['ControlFitnessSummary'])
	if (m<1){
		if (summd>m*summc){type<-"S"} else {type<-"E"}
	} 
	else {
		if (summd>m*summc){type<-"E"} else {type<-"S"}
	}
	return(type)
}

####### Grabs posteriors for each variable for an ORF #######
posmake<-function(orf,orfs,varpos){orfn<-match(orf,orfs)
	norfs<-length(orfs)
	return(lapply(varpos,varposget,orfn,norfs))
}

#### Returns posterior if variable repeated for all ORFs ####
varposget<-function(var,orfn,norfs){
	if (!is.null(dim(var))){z<-var[,orfn]} else {z<-NULL}
	return(z)
}

############################### Likelihood Functions ################################

##### Does max. lik. fit for all colonies, given colonyzer.read or rod.read input #####
qfa.fit<-function(d,inocguess,ORF2gene="ORF2GENE.txt",fmt="%Y-%m-%d_%H-%M-%S",minK=0.025,detectThresh=0.0005,globalOpt=FALSE,logTransform=FALSE,fixG=TRUE,AUCLim=5,STP=20,modelFit=TRUE,...){
	# Define optimization bounds based on inocguess #
	lowK<-max(0.9*inocguess,minK); upK<-1.0
	lowr<-0; upr<-25
	# We often fix inoculation density, but maybe users might prefer to infer it from growth curves
	if(fixG) {lowg<-0.9*inocguess; upg<-1.1*inocguess}else{lowg<-0.01*inocguess; upg<-100.0*inocguess}
	if(globalOpt==TRUE) {lowg<-1e-6*inocguess; upg<-1e6*inocguess}
	lowv<-0.1; upv<-10.0
	xybounds<-list(K=c(lowK,upK),r=c(lowr,upr),g=c(lowg,upg),v=c(lowv,upv))

	# Create orf2gene dictionary
	if (ORF2gene!=FALSE){gdict<-orf2gdict(ORF2gene)}
	# Vector of barcodes
	barcodes<-unique(d$Barcode); nbc<-length(barcodes)
	# Get big data frame ready for results
	results<-data.frame()
	# For each barcode, optimize
	bcount<-0
	for (bcode in barcodes){bcount<-bcount+1
		# Say which plate you're optimizing
		print(paste("Optimizing Plate",bcount,"/",nbc,":",bcode))
		# Restrict data to just this barcode
		dbc<-d[d$Barcode==bcode,]
		inoctime=dbc$Inoc.Time[1]
		# Get list of unique colony positions
		positions<-lapply(1:length(dbc[,1]),index2pos,dbc)
		positions<-unique(positions)
		# Fit logistic model to each colony
		bcfit<-t(sapply(positions,colony.fit,dbc,inocguess,xybounds,globalOpt,detectThresh,minK,logTransform,AUCLim,STP,modelFit,...))
		info<-t(sapply(positions,colony.info,dbc))
		rows<-sapply(positions,rcget,"row")
		cols<-sapply(positions,rcget,"col")
		# Bind Data frame of barcode results to overall results
		# s=bcfit[,4]
		results<-rbind(results,data.frame(Barcode=as.character(info[,1]),Row=rows,Col=cols,
		ScreenID=as.character(info[,9]),Treatment=as.character(info[,2]),Medium=as.character(info[,3]),ORF=as.character(info[,4]),
		K=bcfit[,1],r=bcfit[,2],g=bcfit[,3],v=bcfit[,4],obj=bcfit[,5],t0=bcfit[,6],nAUC=bcfit[,7],nSTP=bcfit[,8],d0=bcfit[,9],Screen.Name=as.character(info[,5]),
		Library.Name=as.character(info[,6]),MasterPlate.Number=as.numeric(info[,7]),Timeseries.order=as.numeric(info[,8]),
		Inoc.Time=inoctime,TileX=as.numeric(info[,10]),TileY=as.numeric(info[,11]),XOffset=as.numeric(info[,12]),YOffset=as.numeric(info[,13]),
		Threshold=as.numeric(info[,14]),EdgeLength=as.numeric(info[,15]),EdgePixels=as.numeric(info[,16]),RepQuad=as.numeric(info[,17])))
	} #bcode
	if (ORF2gene!=FALSE){results$Gene<-sapply(as.character(results$ORF),orf2g,gdict)}
	return(results)
}

makeFitness<-function(results,AUCLim=5,dtmax=25){
	# Fitness definitions from Addinall et al. 2011
	results$MDRMDP=mdrmdp(results$K,results$r,results$g,results$v)	
	results$MDP=mdp(results$K,results$r,results$g,results$v)
	results$MDR=mdr(results$K,results$r,results$g,results$v)
	# Set spurious fitnesses to zero
	results$MDR[(results$r>7)&(results$K<0.0275)]=0
	results$MDRMDP[(results$r>7)&(results$K<0.0275)]=0
	# Doubling time (doublings per hour)
	if("t0"%in%rownames(results)){
		results$DT=dtl(results$K,results$r,results$g,results$v,results$t0)*24
	}else{
		results$DT=(1/results$MDR)*24
	}
	# If doubling time is Inf, set to dtmax
	results$DT[abs(results$DT)>dtmax]=dtmax
	# Area under curve
	AUC<-function(dno,tstar,dat) max(0,integrate(Glogist,lower=0,upper=tstar,K=dat$K[dno],r=dat$r[dno],g=dat$g[dno],v=dat$v[dno],subdivisions=1000)$value - tstar*dat$g[dno])
	results$AUC=sapply(1:length(results[,1]),AUC,tstar=AUCLim,dat=results)
	return(results)
}

# Extract row & col from list of position vectors
rcget<-function(posvec,rc) posvec[match(rc,c("row","col"))] 

# Closure returning a loess-based approximating function for a timeseries
loapproxfun=function(t,g,span=0.2){
	lo=loess(g~t,span=span)
	# If loess smoothing with specified span does not work (e.g. too few points), revert to linear interpolation
	if(is.na(lo$fitted[1])){
		loex=approxfun(t,g,method="linear",rule=2)
	}else{
		loex=function(t){
			cdens=predict(lo,t)
			cdens[t<min(lo$x)]=predict(lo,min(lo$x))
			cdens[t>max(lo$x)]=predict(lo,max(lo$x))
			return(cdens)
		}
	}
	return(loex)
}

### Function that does the optimization for one colony ###
colony.fit<-function(position,bcdata,inocguess,xybounds,globalOpt=FALSE,detectThresh=0,minK=0,logTransform=FALSE,AUCLim=5,STP=10,modelFit=TRUE,...){
	# Get row & column to restrict data
	row<-position[1]; col<-position[2]
	do<-bcdata[(bcdata$Row==row)&(bcdata$Col==col),]
	len1=length(do$Growth)
	# Generate numerical AUC
	if(len1>1){
			loapproxfree=loapproxfun(as.numeric(do$Expt.Time),as.numeric(do$Growth),span=0.5)
			loapprox=function(x) pmax(0,loapproxfree(x))
			nAUC=as.numeric(integrate(loapprox,0,AUCLim)$value)
			nSTP=as.numeric(loapprox(STP))
	}else{
			# If there's only one photograph (e.g. single time point 1536 assay)
			nAUC=NA
			nSTP=do$Growth[1]
	}
	
	# Throw away observations below the detectable threshold
	d=do[as.numeric(do$Growth)>=detectThresh,]
	len2=length(d$Growth)
	# If this has left us with too few points, return "dead colony"
	if ((len2/len1<0.25)|(len2<3)) return(c(inocguess,0,inocguess,1,Inf,0,nAUC,nSTP,do$Growth[1]))
	growth=as.numeric(d$Growth)
	tim=as.numeric(d$Expt.Time)
	maxObs=max(growth)
	# First *detectable* observation at time t0
	t0=min(tim)
	if(modelFit){
		#xybounds$K=c(0.9*maxObs,1.1*maxObs)
		if(globalOpt) {
			pars=de.fit(tim,growth,inocguess,xybounds,initPop=TRUE,logTransform=logTransform)
			# Check for high fraction of K at end of modelled experiment
			GEnd=Glogist(pars[1],pars[2],pars[3],pars[4],max(tim))
			if(GEnd/pars[1]<0.75){ # If experiment not quite finished...
				#print("Modelled growth at end of experiment is less than 0.75 of K estimate...")
				# Put tight bounds on K and optimise again
				Kmin=max(0.95*1.5*GEnd,minK); Kmax=max(1.05*1.5*GEnd,minK)
				xybounds$K=c(Kmin,Kmax)
				pars=de.fit(tim,growth,inocguess,xybounds,initPop=TRUE,logTransform=logTransform)
			}
		}else{
			pars=data.fit(tim,growth,inocguess,xybounds,logTransform=logTransform)
			# Check for high fraction of K at end of modelled experiment
			GEnd=Glogist(pars[1],pars[2],pars[3],pars[4],max(tim))
			if(GEnd/pars[1]<0.75){ # If experiment not quite finished...
				#print("Modelled growth at end of experiment is less than 0.75 of K estimate...")
				# Put tight bounds on K and optimise again
				Kmin=max(0.95*1.5*GEnd,minK); Kmax=max(1.05*1.5*GEnd,minK)
				xybounds$K=c(Kmin,Kmax)
				pars=data.fit(tim,growth,inocguess,xybounds,logTransform=logTransform)
			}
		}
		# Check for spurious combination of relatively high r, low K, spend more time optimising...
		if((mdr(pars[1],pars[2],pars[3],pars[4])>1.0)&((pars[1]<0.05)|(max(growth)<0.05)|(tail(growth,1)<0.05))){ # Try optimising with sick colony as guess
			#print("Attempting to do global fit alternative sick colony growth curve")
			Kmin=max(0.9*inocguess,minK); Kmax=1.5*max(pars[1],minK); xybounds$K=c(Kmin,Kmax); 
			xybounds$r=c(0,3) # Slow growth
			xybounds$v=c(0.75,1.5) # More logistic growth
			inits=list(K=pars[1],r=0.6,g=inocguess,v=1)
			newpars=de.fit(tim,growth,inocguess,xybounds,inits=inits,initPop=TRUE,widenr=FALSE,logTransform=logTransform)			
			if(newpars[5]<=pars[5]) pars=newpars				
		}
		opt=sumsq(pars[1],pars[2],pars[3],pars[4],do$Growth,do$Expt.Time,logTransform=logTransform) # Use all data (not just data below detection thresh)
		dead=sumsq(inocguess,0,inocguess,1,do$Growth,do$Expt.Time,logTransform=logTransform)
		if(dead<=opt) pars=c(inocguess,0,inocguess,1,dead) # Try dead colony
	}else{pars=c(NA,NA,NA,NA,NA)}
	# Add on time of first obs. and numerical AUC, STP
	pars=c(pars,t0,nAUC,nSTP,do$Growth[1])
	return(pars)
}

### Function that fits model to a timecourse
de.fit<-function(tim,growth,inocguess,xybounds,inits=list(),initPop=FALSE,widenr=TRUE,logTransform=FALSE,mxit=2000){
	# Fit to growth curve with differential evolution
	# Get initial guess for parameters
	if(length(inits)==0){
		# Get initial guess for parameters
		init<-guess(tim,growth,inocguess,xybounds)
	}else{init=inits}
	if (widenr) xybounds$r=c(0.01*init$r,10.0*init$r)

	# Function to be optimized
	objf<-function(modpars){
		K<-modpars[1]; r<-modpars[2]
		g<-modpars[3]; v<-abs(modpars[4])
		res=sumsq(K,r,g,v,growth,tim,logTransform=logTransform)
		return(res)
	}

	# Make results repeatable
	#set.seed(1234)
	# Format bounds
	low=c(xybounds$K[1],xybounds$r[1],xybounds$g[1],xybounds$v[1])
	up=c(xybounds$K[2],xybounds$r[2],xybounds$g[2],xybounds$v[2])

	NumParticles=11*length(low)

	if(initPop){
	# Randomly sample from within bounds for initial population
	Klist=runif(NumParticles,min=low[1],max=up[1])
	rlist=runif(NumParticles,min=low[2],max=up[2])
	glist=runif(NumParticles,min=low[3],max=up[3])
	vlist=runif(NumParticles,min=low[4],max=up[4])
	pop=data.frame(K=Klist,r=rlist,g=glist,v=vlist)
	pop=as.matrix(pop)
	# Set first particle equal to initial guess
	pop[1,]=as.numeric(t(init))
	# Set second particle equal to a dead colony
	pop[2,]=c(0.025,0,low[3],1)
	}else{pop=NULL}

	optsol=DEoptim(objf,lower=low,upper=up,
	DEoptim.control(trace = FALSE,NP=NumParticles,initialpop=pop,itermax=mxit)
	)
	pars=abs(as.numeric(optsol$optim$bestmem))
	objval=objf(pars)
	# Sanity check for fitted parameters (no negative growth)
	if (pars[1]<pars[3]){
		pars[1]=pars[3]; pars[2]=0; pars[4]=1
		objval=objf(pars)}
	pars=c(pars,objval)
	return(pars)
}

### Function that fits to for a timecourse
data.fit<-function(tim,growth,inocguess,xybounds,inits=list(),logTransform=FALSE){
	if(length(inits)==0){
		# Get initial guess for parameters
		init<-guess(tim,growth,inocguess,xybounds)
	}else{init=inits}
	#xybounds$r=c(0.75*init$r,1.25*init$r)
	#xybounds$r=c(0.1*init$r,10.0*init$r)
	# Set initial guess to twice best estimate to bias against small r
	#init$r=2.0*init$r
	# Try to stop L-BFGS-B errors by moving away from minK
	xybounds$K[1]=1.01*xybounds$K[1]
	# Function to be optimized
	objf<-function(modpars){
		K<-modpars[1]; r<-modpars[2]
		g<-modpars[3]; v<-abs(modpars[4])
		return(sumsq(K,r,g,v,growth,tim,logTransform=logTransform))
	}
	# Perform optimization
	optsol<-optim(par=unlist(init),fn=objf,gr=NULL,method="L-BFGS-B",
	lower=c(xybounds$K[1],xybounds$r[1],xybounds$g[1],xybounds$v[1]),
	upper=c(xybounds$K[2],xybounds$r[2],xybounds$g[2],xybounds$v[2]),
	control=list(maxit=1000,factr=1e7,trace=0,parscale=c(0.2,10,inocguess,1)))
	pars=abs(as.numeric(optsol$par))
	objval=objf(pars)
	# Sanity check for fitted parameters (no negative growth)
	if (pars[1]<pars[3]){
		pars[1]=pars[3]; pars[2]=0; pars[4]=1
		objval=objf(pars)}
	pars=c(pars,objval)
	if(optsol$message!="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH") print(optsol$message)
	return(pars)
}

### Provide initial guess for logistic model parameters ###
guess<-function(tim,growth,inocguess,xybounds,minK=0.025){
	# Sort time and growth
	growth=growth[order(tim)]
	tim=tim[order(tim)]
	n=length(tim)

	# Enforce positivity and monotonic increasing behaviour in growth
	#growth[1]=max(c(growth[1],0.000000001))
	#for (x in 2:length(growth)) growth[x]=max(c(max(growth[1:(x-1)]),growth[x],0.00000001))
	G0g<-inocguess
	Kg<-max(max(growth),minK)
	vg=1 # Assume logistic model is adequate
	rg=0
	if(n>3){
		# Guess for r is a bit more complicated
		targ<-min(growth)+(max(growth)-min(growth))/2.0
		approxspline=smooth.spline(tim,growth)
		approxcurve=approxfun(approxspline$x,approxspline$y,rule=2)
		solvefn<-function(x) approxcurve(x)-targ
		if((approxcurve(max(tim))>=approxcurve(min(tim)))&(sign(solvefn(min(tim)))!=sign(solvefn(tim[which.max(growth)])))){
			sol=uniroot(solvefn,c(min(tim),tim[which.max(growth)]))
			tmrate=sol$root
		}else{ # Too few points to fit spline curve, still do something sensible...
			tmrate=(min(tim)+max(tim))/2.0;
		}
		if(Kg>G0g) rg=log((Kg-G0g)/G0g)/tmrate
	}
	# Sanity check for guessed parameter values
	# If the data have low correlation, then set r=0
	# If all elements of growth are equal (e.g. zero) then correlation function throws error...
	if((length(unique(growth))==1)|(cor(tim,growth)<0.1)){ rg=0; Kg=minK}
	return(list(K=Kg,r=rg,g=G0g,v=vg))
}#guess

guessNEW<-function(tim,growth,inocguess,xybounds,minK=0.025){
	# Sort time and growth
	growth=growth[order(tim)]
	tim=tim[order(tim)]
	# Enforce positivity and monotonic increasing behaviour in growth
	#growth[1]=max(c(growth[1],0.000000001))
	#for (x in 2:length(growth)) growth[x]=max(c(max(growth[1:(x-1)]),growth[x],0.00000001))
	G0g<-inocguess
	Kg<-max(max(growth),minK)
	vg=1 # Assume logistic model is adequate
	rg=0
	# Subset of data which will be linear on the log scale
	growthadj=growth[growth<0.75*Kg]
	timadj=tim[growth<0.75*Kg]

	n=length(timadj)
	if(n>=1){
		fixed=lm(I(log(growthadj)-log(G0g))~timadj+0)
		fixslope=fixed$coefficients[['timadj']]
		free=lm(log(growthadj)~timadj)
		freeslope=free$coefficients[['timadj']]
		freeint=free$coefficients[['(Intercept)']]
		
		rguess=(freeslope*(Kg/G0g)^vg)/((Kg/G0g)^vg-1)
		if(Kg>G0g) {rg=rguess}else{rg=0}
	}
	# Sanity check for guessed parameter values
	# If the data have low correlation, then set r=0
	# If all elements of growth are equal (e.g. zero) then correlation function throws error...
	if((length(unique(growth))==1)|(cor(tim,growth)<0.1)){ rg=0; Kg=minK}
	return(list(K=Kg,r=rg,g=G0g,v=vg))
}#guess

# Sum of squared error
sumsq<-function(K,r,g,v,growth,tim,logTransform=FALSE){
	if(logTransform){
		ss=sum(((growth-Glogist(K,r,g,v,tim))/growth)^2)/length(growth)
	}else{
		ss=sum((growth-Glogist(K,r,g,v,tim))^2)/length(growth)
	}
	if(is.na(ss)){print("Problem with squared error!"); return(Inf)}else{return(ss)}
}

# Prevent zero or negative growth from occurring
nozero<-function(growth){if (growth<=0){growth<-0.0000000000001} else {growth<-growth}}

#### Get colony information ####
colony.info<-function(position,bcdata){
	# Get row & column to restrict data
	row<-position[1]; col<-position[2]
	d<-bcdata[(bcdata$Row==row)&(bcdata$Col==col),]
	# Want to use the LAST datapoint 
	# Most relevant for edge length (possibly most reliable for position)
	d<-d[order(d$Date.Time,decreasing=TRUE),]
	d<-d[1,]
	# Character vector of colony info
	cvec=c(as.character(d$Barcode),
	as.character(d$Treatments),as.character(d$Medium),as.character(d$ORF),
	as.character(d$Screen.Name),as.character(d$Library.Name),as.character(d$MasterPlate.Number),
	as.character(d$Timeseries.order),as.character(d$ScreenID),as.numeric(d$Tile.Dimensions.X),as.numeric(d$Tile.Dimensions.Y),
	as.numeric(d$X.Offset),as.numeric(d$Y.Offset),as.numeric(d$Threshold),as.numeric(d$Edge.length),as.numeric(d$Edge.Pixels),as.numeric(d$RepQuad))
	if("Client"%in%colnames(d)){cvec$Client=as.character(d$Client)}
	if("ExptDate"%in%colnames(d)){cvec$ExptDate=as.character(d$ExptDate)}
	if("User"%in%colnames(d)){cvec$User=as.character(d$User)}	
	return(cvec)
}

##### Make PDFs #####
qfa.plot<-function(file,results,d,fmt="%Y-%m-%d_%H-%M-%S",barcodes=c(),master.plates=c(),treatments=c(),screen.names=c(),screenIDs=c(),maxg=0,maxt=0,logify=FALSE,densityCol="Growth",curves=TRUE,ylabel="Cell density (AU)",ptype="p"){
	# Sort the data to be plotted sensibly, allowing easy comparison between repeats
	results=results[order(results$MasterPlate.Number,results$Treatment,results$Screen.Name),]
	# Get character vectors of requested barcodes, treatements,etc.; all if none specified
	if (length(master.plates)==0){master.plates<-unique(results$MasterPlate.Number)} 
	results<-results[results$MasterPlate.Number%in%master.plates,]
	d<-d[d$MasterPlate.Number%in%master.plates,]
	if (length(treatments)==0){treatments<-unique(results$Treatment)} 
	results<-results[results$Treatment%in%treatments,]
	d<-d[d$Treatments%in%treatments,]
	if (length(screen.names)==0){screen.names<-unique(results$Screen.Name)} 
	results<-results[results$Screen.Name%in%screen.names,]
	d<-d[d$Screen.Name%in%screen.names,]
	if (length(screenIDs)==0){screenIDs<-unique(results$ScreenID)} 
	results<-results[results$ScreenID%in%screenIDs,]
	d<-d[d$ScreenID%in%screenIDs,]
	if (length(barcodes)==0){barcodes<-unique(results$Barcode)}
	results<-results[results$Barcode%in%barcodes,]
	d<-d[d$Barcode%in%barcodes,]

	print(paste("Plotting for",length(results[,1]),"colonies"))
	# Produce PDF
	colmin<-min(results$Col); rowmin<-min(results$Col)
	colmax<-max(results$Col); rowmax<-max(results$Row)
	pdf(file,4*(colmax-colmin+1),4*(rowmax-rowmin+1))
	if(rowmax*colmax==96) {cexfctr=1.0;marge=c(2,1,2,0.75)} 
	if(rowmax*colmax>96) {cexfctr<-(rowmax*colmax)/384; marge=c(2,1,2,0.75)}
	if(rowmax*colmax>384) {cexfctr<-(rowmax*colmax)/1536; marge=c(2.0,1.5,2.0,1.5)}
	#cexfctr=1.0
	bcount<-0; nbc<-length(barcodes)
	for (bcode in barcodes){
		bcode<-as.character(bcode)
		# Say what you're doing
		bcount<-bcount+1
		print(paste("Plotting for Plate",bcount,"/",nbc,":",bcode))
		# Restricts results and data to this plate
		rbc<-results[as.character(results$Barcode)==bcode,]
		dbc<-d[as.character(d$Barcode)==bcode,]
		# Get starting time for this plate
		inoctime<-rbc$Inoc.Time[1]
		#as.POSIXlt(as.character(bigd$Date.Time),format=fmt)
		# Find number of rows and columns on this plate
		#nrow<-max(rbc$Row)-min(rbc$Row)+1; ncol<-max(rbc$Col)-min(rbc$Col)+1
		nrow<-rowmax-rowmin+1; ncol=colmax-colmin+1
		# Find max growth to set y-axis for this plate
		if (maxg==0) maxg=max(dbc$Growth)
		# Set graphics parameters for each plate
		op<-par(mfrow=c(nrow,ncol),oma=c(13,15,22,1),
		mar=marge,mgp=c(3,1,0),cex=cexfctr)
		## Plot for each row of results for that bcode ##
		z<-apply(rbc,1,rowplot,dbc,inoctime,maxg,fmt,maxt,logify,densityCol=densityCol,curves=curves,ptype=ptype)
		# Title for the plate
		maintit<-paste(rbc$Barcode[1],"Treatment:",rbc$Treatment[1],
		"Medium:",rbc$Medium[1],"Plate:",rbc$MasterPlate.Number[1],sep=" ")
		#cextit<-1008/length(strsplit(maintit,split="")[[1]])
		cextit<-500/nchar(maintit)
		title(main=maintit,xlab="Time since inoculation (days)",line=7,
		ylab=ylabel,cex.main=cextit,cex.lab=8,outer=TRUE)
		  par(op)} #bcode
	dev.off()
}

#### Plot a colony's timecourse from a row of the results #####	
rowplot<-function(resrow,dbc,inoctime,maxg,fmt,maxt,logify,densityCol="Growth",curves=TRUE,ptype="p"){
	row<-as.numeric(resrow['Row']); col<-as.numeric(resrow['Col'])
	if ('Gene'%in%names(resrow)){gene<-resrow['Gene']} else {gene<-resrow['ORF']}
	# Get data for that colony
	dcol<-dbc[(dbc$Row==row)&(dbc$Col==col),]
	growth<-sapply(dcol[[densityCol]],nozero)
	tim<-dcol$Expt.Time
	# Draw the curves and data
	logdraw(row,col,resrow,tim,growth,gene,maxg,maxt=maxt,logify=logify,densityCol=densityCol,curves=curves,ptype=ptype)
}

### Converts row no. to position vector ###	
index2pos<-function(index,dbc) c(dbc[index,'Row'],dbc[index,'Col'])

### Do individual timecourse plot given parameters & data ###
logdraw<-function(row,col,resrow,tim,growth,gene,maxg,fitfunct,maxt=0,scaleT=1.0,logify=FALSE,densityCol="Growth",curves=TRUE,ptype="p"){
	if(logify) {ylog="y"}else{ylog=""}
	plot(NULL,type="n",xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg),log=ylog,xlab="",ylab="",main=gene,frame.plot=0,cex.main=3*scaleT,cex.axis=1*scaleT)
	if(curves){
		# Get logistic parameters, gene name and position
		K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); g<-as.numeric(resrow['g']); v<-as.numeric(resrow['v']);
		MDR<-as.numeric(resrow['MDR']); MDP<-as.numeric(resrow['MDP']); AUC<-as.numeric(resrow['AUC']); DT<-as.numeric(resrow['DT']);
		if(logify) {ylog="y"}else{ylog=""}
		# Add logistic curve
		if(maxt==0) maxt=ceiling(max(tim))
		x=0
		curve(Glogist(K,r,g,v,x),n=31,lwd=2.5,add=TRUE,from=0,to=maxt,xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg))
	}
	# Add data points
	points(tim,growth,col="red",cex=2*scaleT,pch=4,lwd=2,type=ptype,xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg))
	if(curves){
		# Add legend
		legt1<-paste(c("K=","r=","g=","v=","MDR=","MDP=","AUC=","DT="),c(signif(K,3),signif(r,3),signif(g,3),signif(v,3),signif(MDR,3),signif(MDP,3),signif(AUC,3),signif(DT,3)),sep="")
		if(logify){legend("bottomright",legt1,box.lty=0,cex=scaleT)}else{legend("topleft",legt1,box.lty=0,cex=scaleT)}
	}
	legend("topright",sprintf("R%02dC%02d",row,col),box.lty=0,cex=0.5*scaleT)
}

