# -- Library of R functions for Quantile-Frequency Analysis (QFA) --
# by Ta-Hsin Li  (thl@us.ibm.com)  October 30, 2022


#' Trigonometric Quantile Regression (TQR)
#'
#' This function computes trigonometric quantile regression (TQR) for univariate time series at a single frequency.
#' @param y vector of time series
#' @param f0 frequency in [0,1)
#' @param tau sequence of quantile levels in (0,1)
#' @param prepared if \code{TRUE}, intercept is removed and coef of cosine is doubled when \code{f0 = 0.5}
#' @return object of \code{rq()} (coefficients in \code{$coef})
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' fit <- tqr.fit(y,f0=0.1,tau=tau)
#' plot(tau,fit$coef[1,],type='o',pch=0.75,xlab='QUANTILE LEVEL',ylab='TQR COEF')
tqr.fit <- function(y,f0,tau,prepared=TRUE) {

  fix.tqr.coef <- function(coef) {
  # prepare coef from tqr for qdft
  # input:  coef = p*ntau tqr coefficient matrix from tqr.fit()
  # output: 2*ntau matrix of tqr coefficients
    ntau <- ncol(coef)
    if(nrow(coef)==1) {   
      # for f = 0
      coef <- rbind(rep(0,ntau),rep(0,ntau))
    } else if(nrow(coef)==2) {  
      # for f = 0.5: rescale coef of cosine by 2 so qdft can be defined as usual
      coef <- rbind(2*coef[2,],rep(0,ntau))
    } else {
      # for f in (0,0.5)
      coef <- coef[-1,]
    }
    coef
  }
  n <- length(y)
  # create regression design matrix
  tt <- c(1:n)
  if(f0 != 0.5 & f0 != 0) {
    fit <- suppressWarnings(quantreg::rq(y ~ cos(2*pi*f0*tt)+sin(2*pi*f0*tt),tau=tau))
  }
  if(f0 == 0.5) {
    fit <- suppressWarnings(quantreg::rq(y ~ cos(2*pi*f0*tt),tau=tau))
  }
  if(f0 == 0) {
    fit <- suppressWarnings(quantreg::rq(y ~ 1,tau=tau))
  }
  if(prepared) fit$coefficients <- fix.tqr.coef(fit$coefficients)
  fit
}


#' Quantile Discrete Fourier Transform (QDFT)
#'
#' This function computes quantile discrete Fourier transform (QDFT) for univariate or multivariate time series.
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
#' @return matrix or array of the quantile discrete Fourier transform of \code{y}
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y,tau)
#' # Make a cluster for repeated use
#' n.cores <- 2
#' cl <- parallel::makeCluster(n.cores)
#' parallel::clusterExport(cl, c("tqr.fit"))
#' doParallel::registerDoParallel(cl)
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y.qdft <- qdft(y1,tau,n.cores=n.cores,cl=cl)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y.qdft <- qdft(y2,tau,n.cores=n.cores,cl=cl)
#' parallel::stopCluster(cl)
qdft <- function(y,tau,n.cores=1,cl=NULL) {

  z <- function(x) { x <- matrix(x,nrow=2); x[1,]-1i*x[2,] }

  extend.qdft <- function(y,tau,result,sel.f) {
  # define qdft from tqr coefficients
  # input: y = ns*nc-matrix or ns-vector of time series
  #        tau = ntau-vector of quantile levels
  #        result = list of qdft in (0,0.5]
  #        sel.f = index of Fourier frequencies in (0,0.5)
  # output: out = nc*ns*ntau-array or ns*ntau matrix of qdft

    if(!is.matrix(y)) y <- matrix(y,ncol=1)
    nc <- ncol(y)
    ns <- nrow(y)
    ntau <- length(tau)
  
    out <- array(NA,dim=c(nc,ns,ntau))
    for(k in c(1:nc)) {
      # define QDFT at freq 0 as ns * quantile
      out[k,1,] <- ns * quantile(y[,k],tau)
      # retrieve qdft for freq in (0,0.5]
      tmp <- matrix(unlist(result[[k]]),ncol=ntau,byrow=TRUE)
      # define remaining values by conjate symmetry (excluding f=0.5) 
      tmp2 <- NULL
      for(j in c(1:ntau)) {
        tmp2 <- cbind(tmp2,rev(Conj(tmp[sel.f,j])))
      }
      # assemble & rescale everything by ns/2 so that periodogram = |dft|^2/ns
      out[k,c(2:ns),] <- rbind(tmp,tmp2) * ns/2
    }
    if(nc == 1) out <- matrix(out[1,,],ncol=ntau)
    out
  }

  if(!is.matrix(y)) y <- matrix(y,ncol=1)
  ns <- nrow(y)
  nc <- ncol(y)
  f2 <- c(0:(ns-1))/ns
  # Fourier frequencies in (0,0.5]
  f <- f2[which(f2 > 0 & f2 <= 0.5)]
  sel.f <- which(f < 0.5)
  nf <- length(f)
  ntau <- length(tau)
  keep.cl <- TRUE
  if(n.cores>1 & is.null(cl)) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("tqr.fit"))
    doParallel::registerDoParallel(cl)
    keep.cl <- FALSE
  }
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  # compute qdft for f in (0,0.5]
  result<-list()
  i <- 0
  for(k in c(1:nc)) {
    yy <- y[,k] 
    if(n.cores>1) {
      tmp <- foreach::foreach(i=1:nf) %dopar% {
	    tqr.fit(yy,f[i],tau)$coef
      }
    } else {
      tmp <- foreach::foreach(i=1:nf) %do% {
        tqr.fit(yy,f[i],tau)$coef
      }
    }
    # tmp = a list over freq of 2 x ntau coefficiets 
    out <- lapply(tmp,FUN=function(x) {apply(x,2,z)})
    result[[k]] <- out
  }
  if(n.cores>1 & !keep.cl) {
    parallel::stopCluster(cl)
    cl <-NULL
  }
  # extend qdft to f=0 and f in (0.5,1) 
  out <- extend.qdft(y,tau,result,sel.f)
  return(out) 
}


#' Quantile-Frequency Plot
#'
#' This function creates an image plot of quantile spectrum.
#' @param freq sequence of frequencies in (0,0.5) at which quantile spectrum is evaluated
#' @param tau sequence of quantile levels in (0,1) at which quantile spectrum is evaluated
#' @param qper real-valued matrix of quantile spectrum evaluated on the \code{freq} x \code{tau} grid
#' @param rg.qper \code{zlim} for \code{qper} (default = \code{range(qper)})
#' @param rg.tau  \code{ylim} for \code{tau} (default = \code{range(tau)})
#' @param rg.freq \code{xlim} for \code{freq} (default = \code{c(0,0.5)})
#' @param color colors (default = \code{colorRamps::matlab.like2(1024)})
#' @param ylab label of y-axis (default = \code{"QUANTILE LEVEL"})
#' @param xlab label of x-axis (default = \code{"FREQUENCY"})
#' @param tlab title of plot (default = \code{NULL})
#' @param set.par if \code{TRUE}, \code{par()} is set internally (single image)
#' @param legend.plot if \code{TRUE}, legend plot is added
#' @return no return value
#' @import 'fields'
#' @import 'graphics'
#' @import 'colorRamps'
#' @export
qfa.plot <- function(freq,tau,qper,rg.qper=range(qper),rg.tau=range(tau),rg.freq=c(0,0.5),
  color=colorRamps::matlab.like2(1024),ylab="QUANTILE LEVEL",xlab="FREQUENCY",tlab=NULL,
  set.par=TRUE,legend.plot=TRUE) {

  if(set.par) { 
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    graphics::par(mfrow=c(1,1),pty="m",lab=c(7,7,7),mar=c(4,4,3,6)+0.1,las=1)
  }
	
  tmp<-qper
  tmp[tmp<rg.qper[1]]<-rg.qper[1]
  tmp[tmp>rg.qper[2]]<-rg.qper[2]
  graphics::image(freq,tau,tmp,xlab=xlab,ylab=ylab,xlim=rg.freq,ylim=rg.tau,col=color,zlim=rg.qper,frame=F,axes=F)
  graphics::axis(1,line=0.5)
  graphics::axis(2,line=0.5,at=seq(0,1,0.1))
  if(legend.plot) {
    fields::image.plot(zlim=rg.qper,legend.only=TRUE,nlevel = 32,legend.mar = NULL, legend.shrink = 1,col=color)
  }
  if(!is.null(tlab)) graphics::title(tlab)
}


#' Quantile Periodogram and Cross-Periodogram (QPER)
#'
#' This function computes quantile periodogram/cross-periodogram (QPER) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()}
#' @return matrix or array of quantile periodogram/cross-periodogram
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' qper <- qdft2qper(y.qdft)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' qfa.plot(ff[sel.f],tau,Re(qper[sel.f,]))
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qper <- qdft2qper(y.qdft)
#' qfa.plot(ff[sel.f],tau,Re(qper[1,1,sel.f,]))
#' qfa.plot(ff[sel.f],tau,Re(qper[1,2,sel.f,]))
qdft2qper <- function(y.qdft) {

  if(is.matrix(y.qdft)) {
    y.qdft <- array(y.qdft,dim=c(1,dim(y.qdft)[1],dim(y.qdft)[2]))
  }
  nc <- dim(y.qdft)[1]
  ns  <- dim(y.qdft)[2]
  nf  <- ns
  ntau <- dim(y.qdft)[3]
  qper <- array(NA,dim=c(nc,nc,nf,ntau))
  for(k in c(1:nc)) {
     tmp1 <- matrix(y.qdft[k,,],ncol=ntau)
     qper[k,k,,] <- (Mod(tmp1)^2) / ns
     if(k < nc) {
       for(kk in c((k+1):nc)) {
            tmp2 <- matrix(y.qdft[kk,,],ncol=ntau)
            qper[k,kk,,] <- tmp1 * Conj(tmp2) / ns
            qper[kk,k,,] <- Conj(qper[k,kk,,])
       }
     }
  }
  if(nc==1) qper <- matrix(qper[1,1,,],nrow=nf,ncol=ntau)
  return(qper)
}


#' Quantile Autocovariance Function (QACF)
#'
#' This function computes quantile autocovariance function (QACF) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()} or SQDFT from \code{sqdft()}
#' @param return.qser if \code{TRUE}, return quantile series (QSER) along with QACF
#' @return matrix or array of quantile autocovariance function if \code{return.sqer = FALSE} (default), else a list with the following elements:
#'   \item{qacf}{matirx or array of quantile autocovariance function}
#'   \item{qser}{matrix or array of quantile series}
#' @import 'stats'
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' qacf <- qdft2qacf(y.qdft)
#' plot(c(0:9),qacf[c(1:10),1],type='h',xlab="LAG",ylab="QACF")
#' qser <- qdft2qacf(y.qdft,return.qser=TRUE)$qser
#' plot(qser[,1],type='l',xlab="TIME",ylab="QSER")
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qacf <- qdft2qacf(y.qdft)
#' plot(c(0:9),qacf[1,2,c(1:10),1],type='h',xlab="LAG",ylab="QACF")
qdft2qacf <- function(y.qdft,return.qser=FALSE) {

  if(is.matrix(y.qdft)) {
    y.qdft <- array(y.qdft,dim=c(1,dim(y.qdft)[1],dim(y.qdft)[2]))
  }
  nc <- dim(y.qdft)[1]
  ns  <- dim(y.qdft)[2]
  nf  <- ns
  ntau <- dim(y.qdft)[3]
  yy <- array(NA,dim=dim(y.qdft))
  for(k in c(1:nc)) {
    yy[k,,] <- Re(matrix(apply(matrix(y.qdft[k,,],ncol=ntau),2,stats::fft,inverse=TRUE),ncol=ntau)) / ns
  }
  qacf <- array(NA,dim=c(nc,nc,dim(y.qdft)[-1]))
  if(nc > 1) {
    for(j in c(1:ntau)) {
      # demean = TRUE is needed
      tmp <- stats::acf(t(as.matrix(yy[,,j],ncol=nc)), type = "covariance", lag.max = ns-1, plot = FALSE, demean = TRUE)$acf
      for(k in c(1:nc)) {
        for(kk in c(1:nc)) {
          qacf[k,kk,,j] <- tmp[,k,kk] 
	    }
      }
    }
    rm(tmp)
  } else {
    for(j in c(1:ntau)) {
      # demean = TRUE is needed
      tmp <- stats::acf(c(yy[,,j]), type = "covariance", lag.max = ns-1, plot = FALSE, demean = TRUE)$acf 
      qacf[1,1,,j] <- tmp[,1,1]
    }	 
  }
  if(nc==1) {
    qacf <- matrix(qacf[1,1,,],ncol=ntau)
    yy <- matrix(yy[1,,],ncol=ntau)
  }
  if(return.qser) {
    return(list(qacf=qacf,qser = yy))
  } else {
    return(qacf)
  }
}


#' Quantile  Series (QSER)
#'
#' This function computes quantile series (QSER) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()} or SQDFT from \code{sqdft()}
#' @return matrix or array of quantile series
#' @import 'stats'
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' qser <- qdft2qser(y.qdft)
#' plot(qser[,1],type='l',xlab="TIME",ylab="QSER")
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qser <- qdft2qser(y.qdft)
#' plot(qser[1,,1],type='l',xlab="TIME",ylab="QSER")
qdft2qser <- function(y.qdft) {

  if(is.matrix(y.qdft)) {
    y.qdft <- array(y.qdft,dim=c(1,dim(y.qdft)[1],dim(y.qdft)[2]))
  }
  nc <- dim(y.qdft)[1]
  ns  <- dim(y.qdft)[2]
  nf  <- ns
  ntau <- dim(y.qdft)[3]
  qser <- array(NA,dim=dim(y.qdft))
  for(k in c(1:nc)) {
    qser[k,,] <- Re(matrix(apply(matrix(y.qdft[k,,],ncol=ntau),2,stats::fft,inverse=TRUE),ncol=ntau)) / ns
  }
  if(nc==1) {
    qser <- matrix(qser[1,,],ncol=ntau)
  }
  return(qser)
}


#' Lag-Window Estimator of Quantile Spectrum and Cross-Spectrum (QSPEC)
#'
#' This function computes lag-window (LW) estimate of quantile spectrum/cross-spectrum (QSPEC) from QACF.
#' @param qacf matrix or array of QACF from \code{qdft2qacf()}
#' @param M bandwidth parameter of lag window (default = \code{NULL}: quantile periodogram)
#' @return A list with the following elements:
#'   \item{spec}{matrix or array of LW estimate}
#'   \item{lw}{lag-window sequence}
#' @import 'stats'
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' qacf <- qdft2qacf(y.qdft)
#' qper.lw <- qspec.lw(qacf,M=5)$spec
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' qfa.plot(ff[sel.f],tau,Re(qper.lw[sel.f,]))
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qacf <- qdft2qacf(y.qdft)
#' qper.lw <- qspec.lw(qacf,M=5)$spec
#' qfa.plot(ff[sel.f],tau,Re(qper.lw[1,2,sel.f,]))
qspec.lw <- function(qacf,M=NULL) {
    
  if(is.matrix(qacf)) qacf <- array(qacf,dim=c(1,1,nrow(qacf),ncol(qacf)))
  
  nc <- dim(qacf)[1]
  ns <- dim(qacf)[3]
  ntau <- dim(qacf)[4]
  nf <- ns
  nf2 <- 2*ns
    
  Hanning <- function(n,M) {
    # Hanning window
    lags<-c(0:(n-1))
    tmp<-0.5*(1+cos(pi*lags/M))
    tmp[lags>M]<-0
    tmp<-tmp+c(0,rev(c(tmp[-1])))
    tmp
  }

  if(is.null(M)) {
    lw <- rep(1,nf2)
  } else {
    lw<-Hanning(nf2,M)
  }
  qper.lw <- array(NA,dim=c(nc,nc,nf,ntau))
  for(k in c(1:nc)) {
      for(j in c(1:ntau)) {
	    gam <- c(qacf[k,k,,j],0,rev(qacf[k,k,-1,j]))*lw
  	    qper.lw[k,k,,j] <- Mod(stats::fft(gam,inverse=FALSE))[seq(1,nf2,2)]
	  }
      if(k < nc) {
        for(kk in c((k+1):nc)) {
          for(j in c(1:ntau)) {
		    gam <- c(qacf[k,kk,,j],0,rev(qacf[kk,k,-1,j]))*lw
  	        qper.lw[k,kk,,j] <- stats::fft(gam,inverse=FALSE)[seq(1,nf2,2)]
		  }
          qper.lw[kk,k,,] <- Conj(qper.lw[k,kk,,])
        }
	  }
  }
  
  # return spectrum
  if(nc==1) qper.lw <- matrix(qper.lw[1,1,,],ncol=ntau)
  return(out=list(spec=qper.lw,lw=lw))
}


#' Quantile Coherence Spectrum (QCOH)
#'
#' This function computes quantile coherence spectrum (QCOH) from quantile spectrum and cross-spectrum.
#' @param qspec array of quantile spectrum and cross-spectrum
#' @param k index of first series (default = 1)
#' @param kk index of second series (default = 2)
#' @return matrix of quantile coherence evaluated at Fourier frequencies in (0,0.5)
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qacf <- qdft2qacf(y.qdft)
#' qper.lw <- qspec.lw(qacf,M=5)$spec
#' qcoh <- qper2qcoh(qper.lw,k=1,kk=2)
#' qfa.plot(ff[sel.f],tau,Re(qcoh))
qper2qcoh<-function(qspec,k=1,kk=2) {

  nf <- dim(qspec)[3]
  ff <- c(0:(nf-1))/nf
  sel.f <- which(ff > 0 & ff < 0.5)
  coh <- Mod(qspec[k,kk,sel.f,])^2/(Re(qspec[k,k,sel.f,])*Re(qspec[kk,kk,sel.f,]))
  coh[coh > 1] <- 1
  coh
}



qsmooth <- function(qper,link=c('linear','log'),method=c('gamm','sp'),spar="GCV") {
# smooth a vector of real or complex data
# input:   qper = vector of real or complex data to be smoothed
#          link = in linear or log domain for smoothing
#        method = gamm: spline smoothing using gamm()
#                 sp:   spline smoothing using smooth.spline()
#          spar = smoothing parameter (except for GCV) for smooth.spline
# output:  tmp = smoothed data
# imports: gamm() in 'mgcv'
#          corAR1() in 'nlme' (included in 'mgcv')
  
  smooth.spline.gamm <- function(x,y=NULL) {
  # smoothing splines with correlated data using gamm()
  # input:        x = vector of response (if y=NULL) or independent variable
  #               y = vector of response if x is independent variable
    if(is.null(y)) {
      y <- x
	  x <- c(1:length(y)) / length(y)
    }
    # fitting spline with AR(1) residuals
    fit <- mgcv::gamm(y ~ s(x), data=data.frame(x=x,y=y), correlation=nlme::corAR1())$gam 
    return(list(x=x,y=stats::fitted(fit),fit=fit))
  }

  if(spar == "GCV") spar <- NULL
  eps <- 1e-16
  if(is.complex(qper)) {
  	tmp<-Re(qper)
    if(abs(diff(range(tmp))) > 0) {
      if(method[1] == 'gamm') {
        tmp<-smooth.spline.gamm(tmp)$y
      } else {
        tmp<-stats::smooth.spline(tmp,spar=spar)$y
      }
	}
  	tmp2<-Im(qper)	
    if(abs(diff(range(tmp2))) > 0) {
      if(method[1] == 'gamm') {
        tmp2<-smooth.spline.gamm(tmp2)$y
      } else {
        tmp2<-stats::smooth.spline(tmp2,spar=spar)$y
      }
	}
    tmp <- tmp + 1i*tmp2
  } else {
  	tmp<-qper
    if(abs(diff(range(tmp))) > 0) {
      if(link[1]=='log') {
	    sig<-mean(tmp)
	    tmp[tmp<eps*sig] <- eps*sig
        if(method[1] == 'gamm') {
          tmp<-exp(smooth.spline.gamm(log(tmp))$y)	
		} else {
          tmp<-exp(stats::smooth.spline(log(tmp),spar=spar)$y)
		}
      } else {
        if(method[1] == 'gamm') {
          tmp<-smooth.spline.gamm(tmp)$y		
		} else {
          tmp<-stats::smooth.spline(tmp,spar=spar)$y
		}
  	  }
	}
  }
  tmp
}



#' Quantile Smoothing of Quantile Periodogram or Spectral Estimate
#'
#'  This function computes quantile-smoothed version of quantile periodogram/cross-periodogram (QPER) or other quantile spectral estimate.
#' @param qper matrix or array of quantile periodogram/cross-periodogram or spectral estimate
#' @param method smoothing method: \code{"gamm"} for \code{mgcv::gamm()}, \code{"sp"} for \code{stats::smooth.spline()}
#' @param spar smoothing parameter in \code{smooth.spline()} (default = \code{"GCV"})
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
#' @return matrix or array of quantile-smoothed quantile spectral estimate
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @import 'mgcv'
#' @import 'nlme'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' qacf <- qdft2qacf(y.qdft)
#' qper.lw <- qspec.lw(qacf,M=5)$spec
#' qfa.plot(ff[sel.f],tau,Re(qper.lw[1,1,sel.f,]))
#' qper.lwqs <- qsmooth.qper(qper.lw,method="sp",spar=0.9)
#' qfa.plot(ff[sel.f],tau,Re(qper.lwqs[1,1,sel.f,]))
qsmooth.qper <- function(qper,method=c("gamm","sp"),spar="GCV",n.cores=1,cl=NULL) {
  
  if(is.matrix(qper)) qper <- array(qper,dim=c(1,1,nrow(qper),ncol(qper)))
  nc<-dim(qper)[1]
  nf<-dim(qper)[3]  
  ntau<-dim(qper)[4]  
  eps <- 1e-16
  method <- method[1]
  f2 <- c(0:(nf-1))/nf
  # half of the Fourier frequencies excluding 0 for smoothing
  sel.f1 <- which(f2 >= 0 & f2 <= 0.5)
  # half Fourier frequencies excluding 0 and 0.5 for freq folding
  sel.f2 <- which(f2 > 0 & f2 < 0.5)  
  sel.f3 <- which(f2 > 0.5)
  keep.cl <- TRUE
  if(n.cores>1 & is.null(cl)) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("qsmooth"))
    doParallel::registerDoParallel(cl)
    keep.cl <- FALSE
  }
  spar0 <- spar
  spars <- NULL
  qper.sm <-qper
  for(k in c(1:nc)) {
    spars <- rbind(spars,c(k,k,spar0))
	if(!is.null(cl)) {
	  qper.sm[k,k,sel.f1,] <- t(parallel::parApply(cl,Re(qper[k,k,sel.f1,]),1,qsmooth,link="log",method=method,spar=spar0))
	} else {
	  qper.sm[k,k,sel.f1,] <- t(apply(Re(qper[k,k,sel.f1,]),1,qsmooth,link="log",method=method,spar=spar0))
	}
    for(j in c(1:ntau)) {
	  qper.sm[k,k,sel.f3,j] <- rev(qper.sm[k,k,sel.f2,j])
	}    	
    if(k < nc) {
      for(kk in c((k+1):nc)) {
        spars <- rbind(spars,c(k,kk,spar0))
        if(!is.null(cl)) {
          tmp1 <- t(parallel::parApply(cl,Re(qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
          tmp2 <- t(parallel::parApply(cl,Im(qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
        } else {
          tmp1 <- t(apply(Re(qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
          tmp2 <- t(apply(Im(qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
        }
  	    qper.sm[k,kk,,] <- tmp1 + 1i*tmp2
        qper.sm[kk,k,,] <- Conj(qper.sm[k,kk,,])
	  }
	}
  }
  if(n.cores>1 & !keep.cl) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  if(nc==1) qper.sm <- matrix(qper.sm[1,1,,],ncol=ntau)
  return(qper.sm)
}


#' Quantile Smoothing of Quantile Discrete Fourier Transform
#'
#' This function computes quantile-smoothed version of quantile discrete Fourier transform (QDFT).
#' @param y.qdft matrix or array of QDFT from \code{qdft()}
#' @param method smoothing method: \code{"gamm"} for \code{mgcv::gamm()}, \code{"sp"} for \code{stats::smooth.spline()}
#' @param spar smoothing parameter in \code{smooth.spline()} (default = \code{"GCV"})
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = NULL)
#' @return matrix or array of quantile-smoothed QDFT
#' @import 'stats'
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @import 'mgcv'
#' @import 'nlme'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' y.qdft <- qsmooth.qdft(y.qdft,method="sp",spar=0.9)
#' qacf <- qdft2qacf(y.qdft)
#' qper.qslw <- qspec.lw(qacf,M=5)$spec
#' qfa.plot(ff[sel.f],tau,Re(qper.qslw[1,1,sel.f,]))
qsmooth.qdft<-function(y.qdft,method=c('gamm','sp'),spar="GCV",n.cores=1,cl=NULL) {

  if(is.matrix(y.qdft)) y.qdft <- array(y.qdft,dim=c(1,nrow(y.qdft),ncol(y.qdft)))
  nc<-dim(y.qdft)[1]
  nf<-dim(y.qdft)[2] 
  ntau<-dim(y.qdft)[3]
  method <- method[1]
  f2 <- c(0:(nf-1))/nf
  # half of the Fourier frequencies excluding 0 for smoothing
  sel.f1 <- which(f2 > 0 & f2 <= 0.5)
  # half Fourier frequencies excluding 0 and 0.5 for freq folding
  sel.f2 <- which(f2 > 0 & f2 < 0.5)  
  sel.f3 <- which(f2 > 0.5)
  keep.cl <- TRUE
  if(n.cores>1 & is.null(cl)) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("qsmooth"))
    doParallel::registerDoParallel(cl)
    keep.cl <- FALSE
  }
  spar0 <- spar
  spars <- NULL
  qdft.sm <- y.qdft
  for(k in c(1:nc)) {
    spars <- rbind(spars,c(k,spar0))
	if(!is.null(cl)) {
	  tmp1 <- t(parallel::parApply(cl,Re(y.qdft[k,sel.f1,]),1,qsmooth,method=method,spar=spar0))
	  tmp2 <- t(parallel::parApply(cl,Im(y.qdft[k,sel.f1,]),1,qsmooth,method=method,spar=spar0))
	} else {
	  tmp1 <- t(apply(Re(y.qdft[k,sel.f1,]),1,qsmooth,method=method,spar=spar0))
	  tmp2 <- t(apply(Im(y.qdft[k,sel.f1,]),1,qsmooth,method=method,spar=spar0))
	}
	qdft.sm[k,sel.f1,] <- tmp1 + 1i* tmp2	
    for(j in c(1:ntau)) {
	  qdft.sm[k,sel.f3,j] <- Conj(rev(qdft.sm[k,sel.f2,j]))
	}
  }
  if(n.cores>1 & !keep.cl) {
    parallel::stopCluster(cl)
    cl <- NULL
  } 
  if(nc==1) qdft.sm <- matrix(qdft.sm,ncol=ntau)
  return(qdft.sm)
}  


#' Kullback-Leibler Divergence of Quantile Spectral Estimate
#'
#' This function computes Kullback-Leibler divergence (KLD) of quantile spectral estimate.
#' @param qper matrix or array of quantile spectral estimate from, e.g., \code{qspec.lw()}
#' @param qspec matrix of array of true quantile spectrum/cross-spectrum (same dimension as \code{qper})
#' @param sel.f index of selected frequencies for computation (default = \code{NULL}: all frequencies)
#' @param sel.tau index of selected quantile levels for computation (default = \code{NULL}: all quantile levels)
#' @return real number of Kullback-Leibler divergence
#' @export
qkl.divergence <- function(qper,qspec,sel.f=NULL,sel.tau=NULL) {

  if(is.matrix(qper)) {
      qper <- array(qper,dim=c(1,1,nrow(qper),ncol(qper)))
      qspec <- array(qspec,dim=c(1,1,nrow(qspec),ncol(qspec)))
  }
  nc <- dim(qper)[1]
  if(is.null(sel.f)) sel.f <- c(2:dim(qper)[3])
  if(is.null(sel.tau)) sel.tau <- c(2:dim(qper)[4])
  if(nc > 1) {
    out <- 0
    for(j in c(1:length(sel.tau))) {
      for(i in c(1:length(sel.f))) {
         S <- matrix(qper[,,sel.f[i],sel.tau[j]],ncol=nc)
         S0 <- matrix(qspec[,,sel.f[i],sel.tau[j]],ncol=nc)
         tmp1 <- Mod(sum(diag(S %*% solve(S0))))
         tmp2 <- log(prod(abs(Re(diag(qr(S)$qr)))) / prod(abs(Re(diag(qr(S0)$qr)))))
         out <- out + tmp1 - tmp2 - nc
      }
    }
    out <- out / (length(sel.f) * length(sel.tau))
  } else {
    S <- Re(qper[1,1,sel.f,sel.tau])
    S0 <- Re(qspec[,1,sel.f,sel.tau])
    out <- mean(S/S0 - log(S/S0) - 1)    
  }
  out
}


#' @useDynLib qfa, .registration=TRUE
rq.fit.fnb2 <- function(x, y, zeta0, rhs, beta=0.99995,eps=1e-06) {
# modified version of rq.fit.fnb() in 'quantreg' to accept user-provided zeta0 and rhs
# for the Primal-Dual Interior Point algorithm (rqfnb) of Portnoy and Koenker (1997) 
# input:   x = n*p design matrix
#          y = n-vector of response
#      zeta0 = n-vector of initial values of dual variables (x in rqfnb) [ zeta0 = 1-tau in rq.fit.fnb ]
#        rhs = right hand side in the dual formulation [rhs = (1 - tau) * apply(x, 2, sum) in rq.fit.fnb ]
# output: coefficients = p-vector of regression coefficients (primal variables)
#             dualvars = n-vector of dual variables
#                  nit = number of iterations
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x)) 
        stop("x and y don't match n")
    d <- rep(1, n)
    u <- rep(1, n)
    wn <- rep(0,9*n)
    wp <- rep(0,(p+3)*p)
    wn[1:n] <- zeta0  # initial value of dual variables 
    z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d), 
        u = as.double(u), beta = as.double(beta), eps = as.double(eps), 
        wn = as.double(wn), wp = as.double(wp), 
        nit = integer(3), info = integer(1))
    if (z$info != 0) 
        warning(paste("Error info = ", z$info, "in stepy: possibly singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    list(coefficients = coefficients, dualvars = z$wn[(2*n+1):(3*n)], nit = z$nit)
}


#' Spline Quantile Regression (SQR)
#'
#' This function computes the spline quantile regression (SQR) solution given response vector and design matrix.
#' It uses the code \code{rqfnb.f} in the "quantreg" package with the permission of Dr. R. Koenker.
#' @param y response vector
#' @param X design matrix (\code{nrow(X) = length(y)})
#' @param tau sequence of quantile levels in (0,1)
#' @param c0 penalty parameter
#' @param d subsampling rate of quantile levels (default = 1)
#' @param weighted if \code{TRUE}, penalty function is weighted (default = \code{FALSE})
#' @param mthreads if \code{TRUE}, multithread BLAS is enabled when available (default = \code{FALSE}, required for parallel computing) 
#' @return A list with the following elements:
#'   \item{coefficients}{matrix of regression coefficients}
#'   \item{nit}{number of iterations}
#' @import 'stats'
#' @import 'splines'
#' @import 'RhpcBLASctl'
#' @export
sqr.fit <- function(y,X,tau,c0,d=1,weighted=FALSE,mthreads=FALSE) {  

  create_coarse_grid <- function(tau,d=1) {
  # create index of a coarse quantile grid for SQR
    ntau <- length(tau)
    sel.tau0 <- seq(1,ntau,d)
    if(sel.tau0[length(sel.tau0)] < ntau) sel.tau0 <- c(sel.tau0,ntau)
    sel.tau0
  }

  create_weights <- function(tau) {
  # quantile-dependent weights of second derivatives in SQR penalty
    0.25 / (tau*(1-tau))
  }

  create_spline_matrix <- function(tau0,tau) {
  # create spline matrix and its second derivative for SQR
  # input:    tau0 = subset of tau
  #            tau = complete set of quantiles for interpolation
  # output: bsMat  = matrix of basis functions
  #         dbsMat = matrix of second derivatives
  #         bsMat2 = maxtrix of basis function for interpolation to tau
  # imports: 'splines'
    knots <- stats::smooth.spline(tau0,tau0)$fit$knot 
    # rescale tau0 and tau to [0,1] for spline matrix to be consistent with smooth.spline
    bsMat <- splines::splineDesign(knots,(tau0-min(tau0))/(max(tau0)-min(tau0)))
    dbsMat <- splines::splineDesign(knots,(tau0-min(tau0))/(max(tau0)-min(tau0)),derivs=2)
    bsMat2 <- splines::splineDesign(knots,(tau-min(tau))/(max(tau)-min(tau)))
    return(list(bsMat=bsMat,dbsMat=dbsMat,bsMat2=bsMat2,knots=knots))
  }

  rescale_penalty <- function(c0,n,tau0,dbsMat,weighted=FALSE) {
  # rescale penalty par by weights and weighted l1 norm of dbsMat
  # input:   c0 = user-specified penalty parameter
  #           n = number of observations
  #        tau0 = L vector of quantile levels
  #      dbsMat = L*K spline matrix of second derivativies
  #    weighted = TRUE if penalty is weighted
  # dependencies: create_weights()
    L <- length(tau0)
    if(weighted) {
      w <- create_weights(tau0)
    } else {
      w <- rep(1,L)  
    }
    c2 <- c0 * w / sum(abs(w * dbsMat))
    c2 <- c2 * n * L
    c2
  }

  create_Phi_matrix <- function(bsVec,p) {
  # create spline basis matrix for LP and dual LP
  # input: bsVec = K-vector of bs functions (K=number of basis functions)
  #            p = number of parameters
  # output:  Phi = p*pK matrix of basis functions
    K <- length(bsVec)
    Phi <- matrix(0,nrow=p,ncol=p*K)
    for(i in c(1:p)) Phi[i,((i-1)*K+1):(i*K)] <- bsVec
    Phi
  }
  
  # turn-off blas multithread to run in parallel
  if(!mthreads) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  n <- length(y)
  ntau <-length(tau)
  p <- ncol(X)
  # create a coarse quantile grid for LP
  sel.tau0 <- create_coarse_grid(tau,d)
  tau0 <- tau[sel.tau0]
  L <- length(tau0)
  # create spline design matrix with knots provided by smooth.spline
  sdm <- create_spline_matrix(tau0,tau)
  K <- ncol(sdm$bsMat)
  # rescale penalty parameter
  c2 <- rescale_penalty(c0,n,tau0,sdm$dbsMat,weighted)
  # set up the observation vector in dual LP
  y2 <- c(rep(y,L),rep(0,p*L))
  # set up the design matrix and rhs vector in dual LP
  X2 <- matrix(0, nrow=(n+p)*L,ncol=p*K)
  rhs <- rep(0,p*K)
  for(ell in c(1:L)) {
    alpha <- tau0[ell]
    cc <- c2[ell]
    # design matrix for spline coefficients
    X0 <- create_Phi_matrix(sdm$bsMat[ell,],p)
    X0 <- X %*% X0
    sel <- ((ell-1)*n+1):(ell*n)
    X2[sel,] <- X0
    # design matrix (already weighted) for penalty function 
    Z0 <- create_Phi_matrix(sdm$dbsMat[ell,],p)
    # constraint matrix
    sel <- (((ell-1)*p+1):(ell*p))+n*L
    X2[sel,] <- 2 * cc * Z0
    rhs <- rhs + (1-alpha) * apply(X0,2,sum) + cc * apply(Z0,2,sum)
  }
  rm(X0,Z0,sel)
  # default initial value of n*L+p*L dual variables
  zeta0 <- c( rep(1-tau0,each=n), rep(0.5, p*L) )
  # compute the primal-dual LP solution
  fit <- rq.fit.fnb2(x=X2, y=y2, zeta0=zeta0, rhs = rhs)
  theta <- fit$coefficients[c(1:(p*K))]
  zeta <- fit$dualvars
  nit <- fit$nit
  # clean up
  rm(X2,y2,zeta0,rhs,fit)
  # map to interpolated regression solution
  theta <- matrix(theta,ncol=1)
  beta <- matrix(0,nrow=p,ncol=ntau)
  for(ell in c(1:ntau)) {
    Phi <- create_Phi_matrix(sdm$bsMat2[ell,],p)
    beta[,ell] <- Phi %*% theta
  }
  return(list(coefficients=beta,nit=nit))
}


#' Trigonometric Spline Quantile Regression (TSQR)
#'
#' This function computes trigonometric spline quantile regression (TSQR) for univariate time series at a single frequency.
#' @param y vector of time series
#' @param f0 frequency in [0,1)
#' @param tau sequence of quantile levels in (0,1)
#' @param c0 penalty parameter
#' @param d subsampling rate of quantile levels (default = 1)
#' @param weighted if \code{TRUE}, penalty function is weighted (default = \code{FALSE})
#' @param prepared if \code{TRUE}, intercept is removed and coef of cosine is doubled when \code{f0 = 0.5}
#' @return object of \code{sqr.fit()} (coefficients in \code{$coef})
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' fit <- tqr.fit(y,f0=0.1,tau=tau)
#' fit.sqr <- tsqr.fit(y,f0=0.1,tau=tau,c0=0.02,d=4)
#' plot(tau,fit$coef[1,],type='p',xlab='QUANTILE LEVEL',ylab='TQR COEF')
#' lines(tau,fit.sqr$coef[1,],type='l')
tsqr.fit <- function(y,f0,tau,c0,d=1,weighted=FALSE,prepared=TRUE) {
  
  create_trig_design_matrix <- function(f0,n) {
  # create design matrix for trignometric regression of time series
  # input: f0 = frequency in [0,1)
  #         n = length of time series
    tt <- c(1:n)
    if(f0 != 0.5 & f0 != 0) {
      X <- cbind(rep(1,n),cos(2*pi*f0*tt),sin(2*pi*f0*tt))
    }
    if(f0 == 0.5) {
      X <- cbind(rep(1,n),cos(2*pi*0.5*tt))
    }
    if(f0 == 0) {
      X <- matrix(rep(1,n),ncol=1)
    }
    X
  }

  fix_tqr_coef <- function(coef) {
  # prepare coef from tqr for qdft
  # input:  coef = p*ntau tqr coefficient matrix from tqr.fit()
  # output: 2*ntau matrix of tqr coefficients
    ntau <- ncol(coef)
    if(nrow(coef)==1) {   
      # for f = 0
      coef <- rbind(rep(0,ntau),rep(0,ntau))
    } else if(nrow(coef)==2) {  
      # for f = 0.5: rescale coef of cosine by 2 so qdft can be defined as usual
      coef <- rbind(2*coef[2,],rep(0,ntau))
    } else {
      # for f in (0,0.5)
      coef <- coef[-1,]
    }
    coef
  }

  n <- length(y)
  X <- create_trig_design_matrix(f0,n)
  fit <- sqr.fit(y,X,tau,c0=c0,d=d,weighted=weighted)
  if(prepared) fit$coefficients <- fix_tqr_coef(fit$coefficients)
  fit
}


#' Spline Quantile Discrete Fourier Transform (SQDFT)
#'
#' This function computes spline quantile discrete Fourier transform (SQDFT) for univariate or multivariate time series.
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param c0 penalty parameter
#' @param d subsampling rate of quantile levels (default = 1)
#' @param weighted if \code{TRUE}, penalty function is weighted (default = \code{FALSE})
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
#' @return matrix or array of the spline quantile discrete Fourier transform of \code{y}
#' @import 'stats'
#' @import 'splines'
#' @import 'RhpcBLASctl'
#' @import 'quantreg'
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.sqdft <- sqdft(y,tau,c0=0.02,d=4)
#' n <- length(y)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' qacf <- qdft2qacf(y.sqdft)
#' qper.sqrlw <- qspec.lw(qacf,M=5)$spec
#' qfa.plot(ff[sel.f],tau,Re(qper.sqrlw[sel.f,]))
sqdft <- function(y,tau,c0=0.02,d=4,weighted=FALSE,n.cores=1,cl=NULL) {

  z <- function(x) { x <- matrix(x,nrow=2); x[1,]-1i*x[2,] }

  extend_qdft <- function(y,tau,result,sel.f) {
  # define qdft from tqr coefficients
  # input: y = ns*nc-matrix or ns-vector of time series
  #        tau = ntau-vector of quantile levels
  #        result = list of qdft in (0,0.5]
  #        sel.f = index of Fourier frequencies in (0,0.5)
  # output: out = nc*ns*ntau-array or ns*ntau matrix of qdft

    if(!is.matrix(y)) y <- matrix(y,ncol=1)
    nc <- ncol(y)
    ns <- nrow(y)
    ntau <- length(tau)
  
    out <- array(NA,dim=c(nc,ns,ntau))
    for(k in c(1:nc)) {
      # define QDFT at freq 0 as ns * quantile
      out[k,1,] <- ns * stats::quantile(y[,k],tau)
      # retrieve qdft for freq in (0,0.5]
      tmp <- matrix(unlist(result[[k]]),ncol=ntau,byrow=TRUE)
      # define remaining values by conjate symmetry (excluding f=0.5) 
      tmp2 <- NULL
      for(j in c(1:ntau)) {
        tmp2 <- cbind(tmp2,rev(Conj(tmp[sel.f,j])))
      }
      # assemble & rescale everything by ns/2 so that periodogram = |dft|^2/ns
      out[k,c(2:ns),] <- rbind(tmp,tmp2) * ns/2
    }
    if(nc == 1) out <- matrix(out[1,,],ncol=ntau)
    out
  }

  if(!is.matrix(y)) y <- matrix(y,ncol=1)
  ns <- nrow(y)
  nc <- ncol(y)
  f2 <- c(0:(ns-1))/ns
  # compute QR at half of the Fourier frequencies, excluding 0
  f <- f2[which(f2 > 0 & f2 <= 0.5)]
  sel.f <- which(f < 0.5)
  nf <- length(f)
  ntau <- length(tau)
  keep.cl <- TRUE
  if(n.cores>1 & is.null(cl)) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("tsqr.fit","sqr.fit","rq.fit.fnb2"))
    doParallel::registerDoParallel(cl)
    keep.cl <- FALSE
  }
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  i <- 0  
  # qdft for f in (0,0.5]
  result<-list()
  for(k in c(1:nc)) {
    yy <- y[,k] 
    if(n.cores>1) {
        tmp <- foreach::foreach(i=1:nf, .packages='quantreg') %dopar% {
	    coef <- tsqr.fit(yy,f[i],tau,c0=c0,d=d,weighted=weighted)$coef
      }
    } else {
      tmp <- foreach::foreach(i=1:nf) %do% {
        coef <- tsqr.fit(yy,f[i],tau,c0=c0,d=d,weighted=weighted)$coef
      }
    }
    # tmp = a list over freq of 2 x ntau coefficiets 
    out <- lapply(tmp,FUN=function(x) {apply(x,2,z)})
    result[[k]] <- out

  }
  if(n.cores>1 & !keep.cl) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  # extend qdft to f=0 and f in (0.5,1) 
  out <- extend_qdft(y,tau,result,sel.f)
  return(out)
}

