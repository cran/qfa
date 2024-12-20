# -- Library of R functions for Quantile-Frequency Analysis (QFA) --
# by Ta-Hsin Li  (thl024@outlook.com)  October 30, 2022; April 8, 2023; August 21, 2023;
# November 26, 2024; December 14, 2024

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
#' @return matrix or array of quantile discrete Fourier transform of \code{y}
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
  result <- list()
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
    tmp <- lapply(tmp,FUN=function(x) {apply(x,2,z)})
    result[[k]] <- tmp
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
#' @param rqper real-valued matrix of quantile spectrum evaluated on the \code{freq} x \code{tau} grid
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
qfa.plot <- function(freq,tau,rqper,rg.qper=range(rqper),rg.tau=range(tau),rg.freq=c(0,0.5),
  color=colorRamps::matlab.like2(1024),ylab="QUANTILE LEVEL",xlab="FREQUENCY",tlab=NULL,
  set.par=TRUE,legend.plot=TRUE) {

  if(set.par) { 
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    graphics::par(mfrow=c(1,1),pty="m",lab=c(7,7,7),mar=c(4,4,3,6)+0.1,las=1)
  }
	
  tmp<-rqper
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


#' Quantile Periodogram (QPER)
#'
#' This function computes quantile periodogram (QPER) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()}
#' @return matrix or array of quantile periodogram
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' y.qper <- qdft2qper(y.qdft)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' qfa.plot(ff[sel.f],tau,Re(y.qper[sel.f,]))
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' y.qper <- qdft2qper(y.qdft)
#' qfa.plot(ff[sel.f],tau,Re(y.qper[1,1,sel.f,]))
#' qfa.plot(ff[sel.f],tau,Re(y.qper[1,2,sel.f,]))
qdft2qper <- function(y.qdft) {

  if(is.matrix(y.qdft)) {
    y.qdft <- array(y.qdft,dim=c(1,dim(y.qdft)[1],dim(y.qdft)[2]))
  }
  nc <- dim(y.qdft)[1]
  ns  <- dim(y.qdft)[2]
  nf  <- ns
  ntau <- dim(y.qdft)[3]
  y.qper <- array(NA,dim=c(nc,nc,nf,ntau))
  for(k in c(1:nc)) {
     tmp1 <- matrix(y.qdft[k,,],ncol=ntau)
     y.qper[k,k,,] <- (Mod(tmp1)^2) / ns
     if(k < nc) {
       for(kk in c((k+1):nc)) {
            tmp2 <- matrix(y.qdft[kk,,],ncol=ntau)
            y.qper[k,kk,,] <- tmp1 * Conj(tmp2) / ns
            y.qper[kk,k,,] <- Conj(y.qper[k,kk,,])
       }
     }
  }
  if(nc==1) y.qper <- matrix(y.qper[1,1,,],nrow=nf,ncol=ntau)
  return(y.qper)
}



#' Quantile Autocovariance Function (QACF)
#'
#' This function computes quantile autocovariance function (QACF) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()}
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
#' y.qacf <- qdft2qacf(y.qdft)
#' plot(c(0:9),y.qacf[c(1:10),1],type='h',xlab="LAG",ylab="QACF")
#' y.qser <- qdft2qacf(y.qdft,return.qser=TRUE)$qser
#' plot(y.qser[,1],type='l',xlab="TIME",ylab="QSER")
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' y.qacf <- qdft2qacf(y.qdft)
#' plot(c(0:9),y.qacf[1,2,c(1:10),1],type='h',xlab="LAG",ylab="QACF")
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
  y.qacf <- array(NA,dim=c(nc,nc,dim(y.qdft)[-1]))
  if(nc > 1) {
    for(j in c(1:ntau)) {
      # demean = TRUE is needed
      tmp <- stats::acf(t(as.matrix(yy[,,j],ncol=nc)), type = "covariance", lag.max = ns-1, plot = FALSE, demean = TRUE)$acf
      for(k in c(1:nc)) {
        for(kk in c(1:nc)) {
          y.qacf[k,kk,,j] <- tmp[,k,kk] 
	    }
      }
    }
    rm(tmp)
  } else {
    for(j in c(1:ntau)) {
      # demean = TRUE is needed
      tmp <- stats::acf(c(yy[,,j]), type = "covariance", lag.max = ns-1, plot = FALSE, demean = TRUE)$acf 
      y.qacf[1,1,,j] <- tmp[,1,1]
    }	 
  }
  if(nc==1) {
    y.qacf <- matrix(y.qacf[1,1,,],ncol=ntau)
    yy <- matrix(yy[1,,],ncol=ntau)
  }
  if(return.qser) {
    return(list(qacf=y.qacf,qser=yy))
  } else {
    return(y.qacf)
  }
}



#' Quantile Series (QSER)
#'
#' This function computes quantile series (QSER) from QDFT.
#' @param y.qdft matrix or array of QDFT from \code{qdft()}
#' @return matrix or array of quantile series
#' @import 'stats'
#' @export
#' @examples
#' # single time series
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qdft <- qdft(y1,tau)
#' y.qser <- qdft2qser(y.qdft)
#' plot(y.qser[,1],type='l',xlab="TIME",ylab="QSER")
#' # multiple time series
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y.qdft <- qdft(cbind(y1,y2),tau)
#' y.qser <- qdft2qser(y.qdft)
#' plot(y.qser[1,,1],type='l',xlab="TIME",ylab="QSER")
qdft2qser <- function(y.qdft) {

  if(is.matrix(y.qdft)) {
    y.qdft <- array(y.qdft,dim=c(1,dim(y.qdft)[1],dim(y.qdft)[2]))
  }
  nc <- dim(y.qdft)[1]
  ns  <- dim(y.qdft)[2]
  nf  <- ns
  ntau <- dim(y.qdft)[3]
  y.qser <- array(NA,dim=dim(y.qdft))
  for(k in c(1:nc)) {
    y.qser[k,,] <- Re(matrix(apply(matrix(y.qdft[k,,],ncol=ntau),2,stats::fft,inverse=TRUE),ncol=ntau)) / ns
  }
  if(nc==1) {
    y.qser <- matrix(y.qser[1,,],ncol=ntau)
  }
  return(y.qser)
}


# Lag-Window (LW) Estimator of Quantile Spectrum
#
# This function computes lag-window (LW) estimate of quantile spectrum from QACF.
# @param y.qacf matrix or array of pre-calculated QACF from \code{qdft2qacf()}
# @param M bandwidth parameter of lag window (default = \code{NULL}: quantile periodogram)
# @return A list with the following elements:
#   \item{spec}{matrix or array of quantile spectrum}
#   \item{lw}{lag-window sequence}
# @import 'stats'
qacf2speclw <- function(y.qacf,M=NULL) {
  
  if(is.matrix(y.qacf)) y.qacf <- array(y.qacf,dim=c(1,1,dim(y.qacf)))
  
  nc <- dim(y.qacf)[1]
  ns <- dim(y.qacf)[3]
  ntau <- dim(y.qacf)[4]
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
	    gam <- c(y.qacf[k,k,,j],0,rev(y.qacf[k,k,-1,j]))*lw
  	    qper.lw[k,k,,j] <- Mod(stats::fft(gam,inverse=FALSE))[seq(1,nf2,2)]
	  }
      if(k < nc) {
        for(kk in c((k+1):nc)) {
          for(j in c(1:ntau)) {
		    gam <- c(y.qacf[k,kk,,j],0,rev(y.qacf[kk,k,-1,j]))*lw
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




#' Quantile Coherence Spectrum
#'
#' This function computes quantile coherence spectrum (QCOH) from quantile spectrum of multiple time series.
#' @param qspec array of quantile spectrum 
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
#' y.qacf <- qacf(cbind(y1,y2),tau)
#' y.qper.lw <- qspec.lw(y.qacf=y.qacf,M=5)$spec
#' y.qcoh <- qspec2qcoh(y.qper.lw,k=1,kk=2)
#' qfa.plot(ff[sel.f],tau,y.qcoh)
qspec2qcoh<-function(qspec,k=1,kk=2) {

  nf <- dim(qspec)[3]
  ff <- c(0:(nf-1))/nf
  sel.f <- which(ff > 0 & ff < 0.5)
  coh <- Mod(qspec[k,kk,sel.f,])^2/(Re(qspec[k,k,sel.f,])*Re(qspec[kk,kk,sel.f,]))
  coh[coh > 1] <- 1
  coh
}



qsmooth <- function(x,link=c('linear','log'),method=c('gamm','sp'),spar="GCV") {
# smooth a vector of real or complex data
# input:   x = vector of real or complex data to be smoothed
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
  if(is.complex(x)) {
  	tmp<-Re(x)
    if(abs(diff(range(tmp))) > 0) {
      if(method[1] == 'gamm') {
        tmp<-smooth.spline.gamm(tmp)$y
      } else {
        tmp<-stats::smooth.spline(tmp,spar=spar)$y
      }
	}
  	tmp2<-Im(x)	
    if(abs(diff(range(tmp2))) > 0) {
      if(method[1] == 'gamm') {
        tmp2<-smooth.spline.gamm(tmp2)$y
      } else {
        tmp2<-stats::smooth.spline(tmp2,spar=spar)$y
      }
	}
    tmp <- tmp + 1i*tmp2
  } else {
  	tmp<-x
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



# Quantile Smoothing of Quantile Periodogram or Spectral Estimate
#
# This function computes quantile-smoothed version of quantile periodogram or other quantile spectral estimate.
# @param y.qper matrix or array of quantile periodogram or spectral estimate
# @param method smoothing method: \code{"gamm"} for \code{mgcv::gamm()} (default), 
# \code{"sp"} for \code{stats::smooth.spline()}
# @param spar smoothing parameter in \code{smooth.spline()} if \code{method = "sp"} (default = \code{"GCV"})
# @param n.cores number of cores for parallel computing (default = 1)
# @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
# @return matrix or array of quantile-smoothed quantile spectral estimate
# @import 'foreach'
# @import 'parallel'
# @import 'doParallel'
# @import 'mgcv'
# @import 'nlme'
# @export
# @examples
# y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
# y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
# tau <- seq(0.1,0.9,0.05)
# n <- length(y1)
# ff <- c(0:(n-1))/n
# sel.f <- which(ff > 0 & ff < 0.5)
# y.qdft <- qdft(cbind(y1,y2),tau)
# y.qacf <- qdft2qacf(y.qdft)
# y.qper.lw <- qspec.lw(y.qacf=y.qacf,M=5)$spec
# qfa.plot(ff[sel.f],tau,Re(y.qper.lw[1,1,sel.f,]))
# y.qper.lwqs <- qsmooth.qper(y.qper.lw,method="sp",spar=0.9)
# qfa.plot(ff[sel.f],tau,Re(y.qper.lwqs[1,1,sel.f,]))
qsmooth.qper <- function(y.qper,method=c("gamm","sp"),spar="GCV",n.cores=1,cl=NULL) {
  
  if(is.matrix(y.qper)) y.qper <- array(y.qper,dim=c(1,1,dim(y.qper)))
  nc<-dim(y.qper)[1]
  nf<-dim(y.qper)[3]  
  ntau<-dim(y.qper)[4]  
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
  qper.sm <-y.qper
  for(k in c(1:nc)) {
    spars <- rbind(spars,c(k,k,spar0))
	if(!is.null(cl)) {
	  qper.sm[k,k,sel.f1,] <- t(parallel::parApply(cl,Re(y.qper[k,k,sel.f1,]),1,qsmooth,link="log",method=method,spar=spar0))
	} else {
	  qper.sm[k,k,sel.f1,] <- t(apply(Re(y.qper[k,k,sel.f1,]),1,qsmooth,link="log",method=method,spar=spar0))
	}
    for(j in c(1:ntau)) {
	  qper.sm[k,k,sel.f3,j] <- rev(qper.sm[k,k,sel.f2,j])
	}    	
    if(k < nc) {
      for(kk in c((k+1):nc)) {
        spars <- rbind(spars,c(k,kk,spar0))
        if(!is.null(cl)) {
          tmp1 <- t(parallel::parApply(cl,Re(y.qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
          tmp2 <- t(parallel::parApply(cl,Im(y.qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
        } else {
          tmp1 <- t(apply(Re(y.qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
          tmp2 <- t(apply(Im(y.qper[k,kk,,]),1,qsmooth,method=method,spar=spar0))
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



#' Kullback-Leibler Divergence of Quantile Spectral Estimate
#'
#' This function computes Kullback-Leibler divergence (KLD) of quantile spectral estimate.
#' @param y.qper matrix or array of quantile spectral estimate from, e.g., \code{qspec.lw()}
#' @param qspec matrix of array of true quantile spectrum (same dimension as \code{y.qper})
#' @param sel.f index of selected frequencies for computation (default = \code{NULL}: all frequencies)
#' @param sel.tau index of selected quantile levels for computation (default = \code{NULL}: all quantile levels)
#' @return real number of Kullback-Leibler divergence
#' @export
qkl.divergence <- function(y.qper,qspec,sel.f=NULL,sel.tau=NULL) {

  if(is.matrix(y.qper)) {
      y.qper <- array(y.qper,dim=c(1,1,dim(y.qper)))
      qspec <- array(qspec,dim=c(1,1,dim(qspec)))
  }
  nc <- dim(y.qper)[1]
  if(is.null(sel.f)) sel.f <- c(2:dim(y.qper)[3])
  if(is.null(sel.tau)) sel.tau <- c(2:dim(y.qper)[4])
  if(nc > 1) {
    out <- 0
    for(j in c(1:length(sel.tau))) {
      for(i in c(1:length(sel.f))) {
         S <- matrix(y.qper[,,sel.f[i],sel.tau[j]],ncol=nc)
         S0 <- matrix(qspec[,,sel.f[i],sel.tau[j]],ncol=nc)
         tmp1 <- Mod(sum(diag(S %*% solve(S0))))
         tmp2 <- log(prod(abs(Re(diag(qr(S)$qr)))) / prod(abs(Re(diag(qr(S0)$qr)))))
         out <- out + tmp1 - tmp2 - nc
      }
    }
    out <- out / (length(sel.f) * length(sel.tau))
  } else {
    S <- Re(y.qper[1,1,sel.f,sel.tau])
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
#' This function computes spline quantile regression (SQR) solution from response vector and design matrix.
#' It uses the FORTRAN code \code{rqfnb.f} in the "quantreg" package with the kind permission of Dr. R. Koenker.
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




# -- new functions in version 2.0 for SAR and AR estimators of quantile spectra (4/8/2023) --


#' Quantile Autocovariance Function (QACF)
#'
#' This function computes quantile autocovariance function (QACF) from time series 
#' or quantile discrete Fourier transform (QDFT).
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param y.qdft matrix or array of pre-calculated QDFT (default = \code{NULL}: compute from \code{y} and \code{tau});
#' if \code{y.qdft} is supplied, \code{y} and \code{tau} can be left unspecified
#' @param n.cores number of cores for parallel computing of QDFT if \code{y.qdft = NULL} (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing of QDFT (default = \code{NULL})
#' @return matrix or array of quantile autocovariance function
#' @import 'stats'
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' # compute from time series
#' y.qacf <- qacf(y,tau)
#' # compute from QDFT 
#' y.qdft <- qdft(y,tau) 
#' y.qacf <- qacf(y.qdft=y.qdft)
qacf <- function(y,tau,y.qdft=NULL,n.cores=1,cl=NULL) {
  if(is.null(y.qdft)) y.qdft <- qdft(y=y,tau=tau,n.cores=n.cores,cl=cl)
  return(qdft2qacf(y.qdft))
}


#' Quantile Series (QSER)
#'
#' This function computes quantile series (QSER) from time series or quantile discrete Fourier transform (QDFT).
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param y.qdft matrix or array of pre-calculated QDFT (default = \code{NULL}: compute from \code{y} and \code{tau});
#' if \code{y.qdft} is supplied, \code{y} and \code{tau} can be left unspecified
#' @param n.cores number of cores for parallel computing of QDFT if \code{y.qdft = NULL} (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing of QDFT (default = \code{NULL})
#' @return matrix or array of quantile series
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' # compute from time series
#' y.qser <- qser(y,tau)  
#' # compute from QDFT
#' y.qdft <- qdft(y,tau)  
#' y.qser <- qser(y.qdft=y.qdft) 
qser <- function(y,tau,y.qdft=NULL,n.cores=1,cl=NULL) {
  if(is.null(y.qdft)) y.qdft <- qdft(y=y,tau=tau,n.cores=n.cores,cl=cl)
  return(qdft2qser(y.qdft))
}


#' Quantile Periodogram (QPER)
#'
#' This function computes quantile periodogram (QPER) from time series 
#' or quantile discrete Fourier transform (QDFT).
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param y.qdft matrix or array of pre-calculated QDFT (default = \code{NULL}: compute from \code{y} and \code{tau});
#' if \code{y.qdft} is supplied, \code{y} and \code{tau} can be left unspecified
#' @param n.cores number of cores for parallel computing of QDFT if \code{y.qdft = NULL} (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing of QDFT (default = \code{NULL})
#' @return matrix or array of quantile periodogram
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' # compute from time series
#' y.qper <- qper(y,tau)  
#' # compute from QDFT 
#' y.qdft <- qdft(y,tau) 
#' y.qper <- qper(y.qdft=y.qdft)  
qper <- function(y,tau,y.qdft=NULL,n.cores=1,cl=NULL) {
  if(is.null(y.qdft)) y.qdft <- qdft(y=y,tau=tau,n.cores=n.cores,cl=cl)
  return(qdft2qper(y.qdft))
}



# Quantile Spectrum from AR Model of Quantile Series
#
# This function computes quantile spectrum (QSPEC) from an AR model of Quantile Series (QSER).
# @param fit object of AR model from \code{qser2sar()} or \code{qser2ar()}
# @param freq sequence of frequencies in [0,1) (default = \code{NULL}: all Fourier frequencies)
# @return a list with the following elements:
#   \item{spec}{matrix or array of quantile spectrum}
#   \item{freq}{sequence of frequencies}
# @export
ar2qspec<-function(fit,freq=NULL) {
  V <- fit$V
  A <- fit$A
  n <- fit$n
  p <- fit$p
  if(is.vector(V)) {
    V <- array(V,dim=c(1,1,length(V)))
	if(p > 0) A <- array(A,dim=c(1,1,dim(A)))
  }
  nc <- dim(V)[1]
  ntau <- dim(V)[3]
  if(is.null(freq)) freq<-c(0:(n-1))/n
  nf <- length(freq)
  spec <- array(0,dim=c(nc,nc,nf,ntau))
  for(ell in c(1:ntau)) {
    for(jj in c(1:nf)) {
	  if(p > 0) {
        U <- matrix(0,nc,nc)
        for(j in c(1:p)) {
          U <- U + A[,,j,ell] * exp(-1i*2*pi*j*freq[jj])
        }
        U <- matrix(solve(diag(nc) - U),ncol=nc)
        spec[,,jj,ell] = U %*% (matrix(V[,,ell],nrow=nc,ncol=nc) %*% t(Conj(U)))
      } else {
	    spec[,,jj,ell] = matrix(V[,,ell],nrow=nc,ncol=nc)
	  }
    }
  }
  if(nc==1) spec <- spec[1,1,,]
  return(list(spec=spec, freq=freq))
}



#' Spline Autoregression (SAR) Model of Quantile Series
#'
#' This function fits spline autoregression (SAR) model to quantile series (QSER).
#' @param y.qser matrix or array of pre-calculated QSER, e.g., using \code{qser()}
#' @param tau sequence of quantile levels where \code{y.qser} is calculated
#' @param d subsampling rate of quantile levels (default = 1)
#' @param p order of SAR model (default = \code{NULL}: automatically selected by AIC)
#' @param order.max maximum order for AIC if \code{p = NULL} (default = \code{NULL}: determined by \code{stats::ar()})
#' @param spar penalty parameter alla \code{smooth.spline} (default = \code{NULL}: automatically selected)
#' @param method criterion for penalty parameter selection:  \code{"AIC"} (default), \code{"BIC"}, or \code{"GCV"}
#' @param weighted if \code{TRUE}, penalty function is weighted (default = \code{FALSE})
#' @return a list with the following elements:
#'   \item{A}{matrix or array of SAR coefficients}
#'   \item{V}{vector or matrix of SAR residual covariance}
#'   \item{p}{order of SAR model}
#'   \item{spar}{penalty parameter}
#'   \item{tau}{sequence of quantile levels}
#'   \item{n}{length of time series}
#'   \item{d}{subsampling rate of quantile levels}
#'   \item{weighted}{option for weighted penalty function}
#'   \item{fit}{object containing details of SAR fit}
#' @import 'stats'
#' @import 'splines'
#' @export
qser2sar <- function(y.qser,tau,d=1,p=NULL,order.max=NULL,spar=NULL,method=c("GCV","AIC","BIC"),weighted=FALSE) {

  create_coarse_grid <- function(tau,d=1) {
  # create index of a coarse quantile grid for SQR
    ntau <- length(tau)
    sel.tau0 <- seq(1,ntau,d)
    if(sel.tau0[length(sel.tau0)] < ntau) sel.tau0 <- c(sel.tau0,ntau)
    sel.tau0
  }

  create_weights <- function(tau,weighted=TRUE) {
  # quantile-dependent weights for the penalty
    if(weighted) {
      0.25 / (tau*(1-tau))
    } else {
      rep(1,length(tau))
    }
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

  sar_rescale_penalty <- function(Z2,D2) {
  # rescalng factor for the SAR penalty parameter alla smooth.spline
    r <- sum(diag(Z2))
    r <- r / sum(diag(D2))
	r
  }

  sar_obj <- function(spar,Z2,D2,YZ,Y,Z,r,nc,L,n,p,method=c("GCV","AIC","BIC"),return.solution=FALSE) {
  # this function computes the criterion for penalty parameter selection in SAR using spar
  # with the option of returning the details of fit
    lam <- r * 256^(3*spar-1)
    # compute least-squares solution
    Z2 <- Z2 + lam * D2
    Z2 <- solve(Z2)
    Theta <- YZ %*% Z2
    # compute residual covariance matrix
    df <- 0
    rss <- 0
    V0 <- array(0,dim=c(nc,nc,L))
    for(l in c(1:L)) {
      df <- df + sum(diag(crossprod(Z[[l]],Z2) %*% Z[[l]])) * nc
      resid <- Y[[l]] - Theta %*% Z[[l]]
      if(method[1]=="GCV") {
        rss <- rss + sum(resid^2)
      }
	  if(method[1] %in% c("AIC","BIC") | return.solution) {
        V0[,,l] <- matrix(tcrossprod(resid),ncol=nc) / n
      }
    }
    # compute parameter selection criterion
    crit <- 0
    if(method[1]=="AIC") {
      for(l in c(1:L)) {
        if(nc==1) {
          crit <- crit + c(abs(V0[,,l])) 
        } else {
          crit <- crit + c(determinant(V0[,,l])$modulus)
        }
      }
      crit <- n*crit + 2*df
    }
    if(method[1]=="BIC") {
      for(l in c(1:L)) {
        if(nc==1) {
          crit <- crit + c(abs(V0[,,l]))
        } else {
          crit <- crit + c(determinant(V0[,,l])$modulus)
        }
      }
      crit <- n*crit + log(n)*df
    }
	if(method[1]=="GCV") {
      crit <- (rss/(L*(n-p))) / (1-df/(L*(n-p)))^2
    }
    if(return.solution) {
      return(list(crit=crit,Theta=Theta,V0=V0,lam=lam,df=df,method=method[1]))
    } else { 
      return(crit)
    }
  }

  if( is.null(spar) & !(method[1] %in% c("GCV","AIC","BIC")) ) {
    stop("method of smoothing parameter selection not GCV, AIC, or BIC!")
  }

  if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))
  nc <- dim(y.qser)[1]
  n <- dim(y.qser)[2]
  ntau <- dim(y.qser)[3]
  
  if( ntau < 10 ) {
    stop("too few quantile levels (must be 10 or more)!")
  }
  if( ntau/d < 10 ) {
    stop("downsampling rate d too large!")
  }
  if( is.null(p) ) {
    if( is.null(order.max) ) {
       stop("order.max unspecified!")	
	} else {
      if( order.max >= n/2 ) {
        stop("order.max too large (must be less than n/2)!")
      }
    }
  }
  if( !is.null(p) ) {
    if(p >= n) {
      stop("order p too large (must be less than n)!")
    }
  }

  # create a coarse quantile grid for SAR
  sel.tau0 <- create_coarse_grid(tau,d)
  tau0 <- tau[sel.tau0]
  L <- length(tau0)
  # create L-by-K spline design matrix with knots provided by smooth.spline
  sdm <- create_spline_matrix(tau0,tau)
  K <- ncol(sdm$bsMat)
  # create weights
  w <- create_weights(tau0,weighted)
  
  # fit VAR to obtain AIC for order selection
  if(is.null(p)) {
    aic <- NULL
    for(l in c(1:L)) {
      mu <- apply(matrix(y.qser[,,sel.tau0[l]],nrow=nc),1,mean)
      aic <- cbind(aic,stats::ar(t(matrix(y.qser[,,sel.tau0[l]]-mu,nrow=nc)),
                   order.max = order.max, method=c("yule-walker"))$aic)
    }
    aic <- apply(aic,1,mean)
    # optimal order minimizes average aic
    p <- as.numeric(which.min(aic)) - 1
  }

  # SAR model of order 0
  if(p == 0) {
     A <- NULL	
     V0 <- array(0,dim=c(nc,nc,L)) 
     V <- array(0,dim=c(nc,nc,ntau))
     for(ell in c(1:ntau)) {
       mu <- apply(matrix(y.qser[,,ell],nrow=nc),1,mean)
       V[,,ell] <- matrix(tcrossprod(matrix(y.qser[,,ell]-mu,nrow=nc)),ncol=nc) / n
     }
     if(nc==1) {
	   V0 <- V0[1,1,]
       V <- V[1,1,]
     }
     return(list(A=A,V=V,p=p,spar=NULL,tau=tau,n=n,d=d,weighted=weighted,fit=NULL))
  }
  
  # SAR model of order p > 0 
  Z2 <- matrix(0,nrow=K*nc*p,ncol=K*nc*p)
  D2 <- matrix(0,nrow=K*nc*p,ncol=K*nc*p)
  YZ <- matrix(0,nrow=nc,ncol=K*nc*p)
  Z <- list()
  Y <- list()
  for(l in c(1:L)) {
    Zt <- matrix(0,nrow=n-p,ncol=K*nc*p)
    mu <- apply(matrix(y.qser[,,sel.tau0[l]],nrow=1),1,mean)
    for(j in c(1:p)) {
      Yt <- t(matrix(y.qser[,c((p-j+1):(n-j)),sel.tau0[l]]-mu,nrow=nc))
      j0 <- (j-1)*K*nc
      for(k in c(1:K)) {
        Zt[,(j0+(k-1)*nc+1):(j0+k*nc)] <- sdm$bsMat[l,k] * Yt
      }
    }
    D0 <- matrix(0,nrow=nc,ncol=K*nc)
    for(k in c(1:K)) {
      D0[,((k-1)*nc+1):(k*nc)] <- diag(sdm$dbsMat[l,k],nrow=nc,ncol=nc)
    }
	Dt <- matrix(0,nrow=nc*p,ncol=K*nc*p)	
    for(j in c(1:p)) {
      Dt[((j-1)*nc+1):(j*nc),((j-1)*K*nc+1):(j*K*nc)] <- D0
    }
    Yt <- t(matrix(y.qser[,c((p+1):n),sel.tau0[l]]-mu,nrow=nc))
    Z2 <- Z2 + crossprod(Zt)
    YZ <- YZ + crossprod(Yt,Zt)
    D2 <- D2 + w[l] * crossprod(Dt)
    Z[[l]] <- t(Zt)
    Y[[l]] <- t(Yt)
  }
  # clean up
  rm(Zt,Dt,Yt,D0)
  # compute scaling factor of penalty parameter alla smooth.spline
  r <- sar_rescale_penalty(Z2,D2)
  # select optimal penalty parameter
  if(is.null(spar)) {
    spar <- stats::optimize(sar_obj,interval=c(-1.5,1.5),
	        Z2=Z2,D2=D2,YZ=YZ,Y=Y,Z=Z,r=r,nc=nc,L=L,n=n,p=p,method=method[1])$minimum
  }
  # compute SAR solution
  fit <- sar_obj(spar,Z2,D2,YZ,Y,Z,r,nc,L,n,p,method=method[1],return.solution=TRUE)
  # clean up    
  rm(Z,Y,YZ,Z2,D2)
  # compute SAR coefficient functions
  A <- array(0,dim=c(nc,nc,p,ntau))
  for(j in c(1:p)) {
    theta <- matrix(fit$Theta[,((j-1)*K*nc+1):(j*K*nc)],nrow=nc)
    for(ell in c(1:ntau)) {
	  Phi <- matrix(0,nrow=K*nc,ncol=nc)
      for(k in c(1:K)) {
        Phi[((k-1)*nc+1):(k*nc),] <- diag(sdm$bsMat2[ell,k],nrow=nc,ncol=nc)
      }
      A[,,j,ell] <- theta %*% Phi
    }
  }
  # interpolate residual covariance matrice
  V <- array(0,dim=c(nc,nc,ntau))
  for(k in c(1:nc)) {
    fit.ss <- stats::smooth.spline(tau0,fit$V0[k,k,],lambda=fit$lam)
    V[k,k,] <- stats::predict(fit.ss,tau)$y
    if(k < nc) {
      for(kk in c((k+1):nc)) {
        fit.ss <- stats::smooth.spline(tau0,fit$V0[k,kk,],lambda=fit$lam)
        V[k,kk,] <- stats::predict(fit.ss,tau)$y
	  }
      V[kk,k,] <- V[k,kk,]
    }
  }
  # fix possible nonpositive definite matrices
  eps <- 1e-8
  for(ell in c(1:ntau)) {
    while(min(eigen(V[,,ell])$values) < 0.0) {
      V[,,ell] <- V[,,ell] + diag(eps,nrow=nc,ncol=nc)
    }
  }
  # special case of nc=1
  if(nc==1) {
    A <- matrix(A[1,1,,],nrow=p)
    V <- V[1,1,]
    fit$V0 <- fit$V0[1,1,]
  }
  return(list(A=A,V=V,p=p,spar=spar,tau=tau,n=n,d=d,weighted=weighted,fit=fit))
}



#' Spline Autoregression (SAR) Estimator of Quantile Spectrum
#'
#' This function computes spline autoregression (SAR) estimate of quantile spectrum.
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param y.qser matrix or array of pre-calculated QSER (default = \code{NULL}: compute from \code{y} and \code{tau});
#' if \code{y.qser} is supplied, \code{y} can be left unspecified
#' @param tau sequence of quantile levels in (0,1)
#' @param d subsampling rate of quantile levels (default = 1)
#' @param p order of SAR model (default = \code{NULL}: automatically selected by AIC)
#' @param order.max maximum order for AIC if \code{p = NULL} (default = \code{NULL}: determined by \code{stats::ar()})
#' @param spar penalty parameter alla \code{smooth.spline} (default = \code{NULL}: automatically selected)
#' @param method criterion for penalty parameter selection: \code{"GCV"}, \code{"AIC"} (default), or \code{"BIC"}
#' @param weighted if \code{TRUE}, penalty function is weighted (default = \code{FALSE})
#' @param freq sequence of frequencies in [0,1) (default = \code{NULL}: all Fourier frequencies)
#' @param n.cores number of cores for parallel computing of QDFT if \code{y.qser = NULL} (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing of QDFT (default = \code{NULL})
#' @return a list with the following elements:
#'   \item{spec}{matrix or array of SAR quantile spectrum}
#'   \item{freq}{sequence of frequencies}
#'   \item{fit}{object of SAR model}
#'   \item{qser}{matrix or array of quantile series if \code{y.qser = NULL}}
#' @import 'stats'
#' @import 'splines'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' # compute from time series
#' y.sar <- qspec.sar(cbind(y1,y2),tau=tau,p=1)
#' qfa.plot(ff[sel.f],tau,Re(y.sar$spec[1,1,sel.f,]))
#' # compute from quantile series
#' y.qser <- qser(cbind(y1,y2),tau)
#' y.sar <- qspec.sar(y.qser=y.qser,tau=tau,p=1)
#' qfa.plot(ff[sel.f],tau,Re(y.sar$spec[1,1,sel.f,]))
qspec.sar <- function(y,y.qser=NULL,tau,d=1,p=NULL,order.max=NULL,spar=NULL,method=c("GCV","AIC","BIC"),
weighted=FALSE,freq=NULL,n.cores=1,cl=NULL) {
  return.qser <- FALSE
  if(is.null(y.qser)) { 
    y.qser <- qser(y,tau,n.cores=n.cores,cl=cl)
    return.qser <- TRUE
  }
  fit <- qser2sar(y.qser,tau=tau,d=d,p=p,order.max=order.max,spar=spar,method=method[1],weighted=weighted)
  tmp <- ar2qspec(fit)
  if(return.qser) {
    return(list(spec=tmp$spec,freq=tmp$freq,fit=fit,qser=y.qser)) 
  } else {
    return(list(spec=tmp$spec,freq=tmp$freq,fit=fit))
  }
}


#' Quantile Periodogram Type II (QPER2)
#'
#' This function computes type-II quantile periodogram for univariate time series.
#' @param y univariate time series
#' @param freq sequence of frequencies in [0,1)
#' @param tau sequence of quantile levels in (0,1)
#' @param weights sequence of weights in quantile regression (default = \code{NULL}: weights equal to 1)
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
#' @return matrix of quantile periodogram evaluated on \code{freq * tau} grid
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @import 'quantreg'
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' y.qper2 <- qper2(y,ff,tau)
#' qfa.plot(ff[sel.f],tau,Re(y.qper2[sel.f,]))
qper2 <- function(y,freq,tau,weights=NULL,n.cores=1,cl=NULL) {

 rsoid2 <- function(tt,ff,a,b) {
 # this function computes trigonometric regression function
  one <- rep(1,length(tt))
  tmp <- (one %o% a) * cos(2*pi*(tt %o% ff)) + (one %o% b) * sin(2*pi*(tt %o% ff))
  tmp <- apply(tmp,1,sum)
  return(tmp)
 }

 qr.cost <- function(y,tau,weights=NULL) {
 # this function computes the cost of quantile regression
   n <- length(y)
   if(is.null(weights)) weights<-rep(1,n)
   tmp <- tau*y
   sel <- which(y < 0)
   if(length(sel) > 0) tmp[sel] <- (tau-1)*y[sel]
   return(sum(tmp*weights,na.rm=T))
 }
 
 tqr.cost <- function(yy,ff,tau,weights=NULL,method='fn') {
 # this function computes the cost of trigonometric quantile regression
 # for a vector of quantile levels tau at a single frequency ff
    n <- length(yy)
    tt <- c(1:n)
	ntau <- length(tau)
    if(ff == 0.5) {
      fit <- suppressWarnings(quantreg::rq(yy ~ cos(2*pi*ff*tt),method=method,tau=tau,weights=weights))
      fit$coefficients <- rbind(fit$coefficients,rep(0,ntau))
    }
    if(ff == 0) {
	  fit <- suppressWarnings(quantreg::rq(yy ~ 1,method=method,tau=tau,weights=weights))
      fit$coefficients <- rbind(fit$coefficients,rep(0,ntau),rep(0,ntau))
    }
    if(ff != 0.5 & ff != 0) {
      fit <- suppressWarnings(quantreg::rq(yy ~ cos(2*pi*ff*tt)+sin(2*pi*ff*tt),method=method,tau=tau,weights=weights))
    }
    fit$coefficients <- matrix(fit$coefficients,ncol=ntau)
    if(ntau == 1) fit$coefficients <- matrix(fit$coefficients,ncol=1)
    cost <- rep(NA,ntau)
    for(i.tau in c(1:ntau)) {
       resid <- yy - fit$coefficients[1,i.tau]-rsoid2(tt,ff,fit$coefficients[2,i.tau],fit$coefficients[3,i.tau])
       cost[i.tau] <- qr.cost(resid,tau[i.tau],weights)
    }
    return(cost)
 }
 
 n <- length(y)
 nf <- length(freq)
 ntau <- length(tau)
 if(is.null(weights)) weights<-rep(1,n)
 keep.cl <- TRUE
 if(n.cores>1 & is.null(cl)) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("rq"))
    doParallel::registerDoParallel(cl)
    keep.cl <- FALSE
 }
 `%dopar%` <- foreach::`%dopar%`
 `%do%` <- foreach::`%do%`
 # cost without trigonometric regressor
 cost0 <- tqr.cost(y,0,tau,weights)
 # cost with trigonometric regressor
 i <- 0  
 if(n.cores>1) {
   out <- foreach::foreach(i=1:nf) %dopar% { tqr.cost(y,freq[i],tau,weights) }
 } else {
   out <- foreach::foreach(i=1:nf) %do% { tqr.cost(y,freq[i],tau,weights) }
 }
 out <- matrix(unlist(out),ncol=ntau,byrow=T)
 out <- matrix(rep(cost0,nf),ncol=ntau,byrow=T) - out
 if(n.cores>1 & !keep.cl) {
    parallel::stopCluster(cl)
    cl <-NULL
 }
 out[out<0]<-0
 if(ntau==1) out<-c(out)
 return(out)
}

  
#' Extraction of SAR Coefficients for Granger-Causality Analysis
#'
#' This function extracts the spline autoregression (SAR) coefficients from an SAR model for Granger-causality analysis.
#' See \code{sar.gc.bootstrap} for more details regarding the use of \code{index}.
#' @param fit object of SAR model from \code{qser2sar()} or \code{qspec.sar()$fit}
#' @param index a pair of component indices for multiple time series 
#' or a sequence of lags for single time series (default = \code{c(1,2)})
#' @return matrix of selected SAR coefficients (number of lags by number of quantiles)
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.sar <- qspec.sar(cbind(y1,y2),tau=tau,p=1)
#' A <- sar.gc.coef(y.sar$fit,index=c(1,2))
sar.gc.coef <- function(fit,index=c(1,2)) {

  test_for_validity <- function(index,nc,p) {
    if(p==0) {
      stop("method not applicable to order 0 models!")	  
    }
    if(nc==1) {
      if(sum(index < 1 | index > p) > 0) {
        stop(paste0("index outside the range of lags 1-",p,"!"))   
      }
    } else {
      if(length(index) != 2) {
        stop(paste0("length of index not equal 2!"))	  	
      }
      if(sum(index < 1 | index > nc) > 0) {
        stop(paste0("index outside the range of time series components 1-",nc,"!"))	   
      }
    }
  }

  if(is.vector(fit$V)) {
	nc <- 1
  } else {
    nc <- dim(fit$V)[1]	  
  }
  p <- fit$p
  test_for_validity(index,nc,p)
  if(nc == 1) {
    return(matrix(fit$A[index,],nrow=length(index)))
  } else {
    return(matrix(fit$A[index[1],index[2],,],nrow=p))
  }
}

  
#' Bootstrap Simulation of SAR Coefficients for Granger-Causality Analysis
#'
#' This function simulates bootstrap samples of selected spline autoregression (SAR) coefficients 
#' for Granger-causality analysis based on the SAR model of quantile series (QSER) under H0: 
#' (a) for multiple time series, the second series specified in \code{index} is not causal 
#' for the first series specified in \code{index};
#' (b) for single time series, the series is not causal at the lags specified in \code{index}.
#' @param y.qser matrix or array of QSER from \code{qser()} or \code{qspec.sar()$qser}
#' @param fit object of SAR model from \code{qser2sar()} or \code{qspec.sar()$fit}
#' @param index a pair of component indices for multiple time series 
#' or a sequence of lags for single time series (default = \code{c(1,2)})
#' @param nsim number of bootstrap samples (default = 1000)
#' @param method method of residual calculation: \code{"ar"} (default) or \code{"sar"}
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param mthreads if \code{TRUE}, multithread BLAS is enabled when available (default = \code{FALSE}, required for parallel computing) 
#' @param seed seed for random sampling (default = \code{1234567})
#' @return array of simulated bootstrap samples of selected SAR coefficients
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @import 'RhpcBLASctl'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.sar <- qspec.sar(cbind(y1,y2),tau=tau,p=1)
#' A.sim <- sar.gc.bootstrap(y.sar$qser,y.sar$fit,index=c(1,2),nsim=5)
sar.gc.bootstrap <- function(y.qser,fit,index=c(1,2),nsim=1000,method=c("ar","sar"),n.cores=1,mthreads=FALSE,seed=1234567) {

  simulate_qser_array <- function(resid,sample.idx,fit,idx0) {
  # this function simulates quantile series under H0 from bootstrap samples of residuals
  # input: resid = nc*(n-p)*ntau-array of residuals from qser2ar()
  #        sample.idx = nn-vector of sampled time indices
  #        fit = SAR model from qser2sar()  
  #        idx0 = indices of SAR coefficients to set to zero under H0
  # output: nc*n*ntau array of simulated quantile series
      p <- fit$p
      nc <- dim(resid)[1]
	  n <- dim(resid)[2] + p
      ntau <- dim(resid)[3]
      nn <- length(sample.idx)
      resid.sim <- array(resid[,sample.idx,],dim=c(nc,nn,ntau))
      qser.sim <- array(0,dim=c(nc,nn,ntau))
      for(ell in c(1:ntau)) {
        for(jj in c((p+1):nn)) {
          if(nc == 1) {
            for(j in c(1:p)) {
              A0 <- fit$A[j,ell]
              if(j %in% idx0) A0 <- 0
              A0 <- matrix(A0,ncol=nc,nrow=nc)			  
              qser.sim[,jj,ell] <- qser.sim[,jj,ell] + A0 %*% matrix(qser.sim[,jj-j,ell],nrow=nc)
            }		  
          } else {
            for(j in c(1:p)) {
              A0 <- fit$A[,,j,ell]
              A0[idx0[1],idx0[2]] <- 0
              qser.sim[,jj,ell] <- qser.sim[,jj,ell] + A0 %*% matrix(qser.sim[,jj-j,ell],nrow=nc)
            } 
          }
          qser.sim[,jj,ell] <- qser.sim[,jj,ell] + resid.sim[,jj,ell] 
        }
      }
      qser.sim <- qser.sim[,c((nn-n+1):nn),]
      return(qser.sim)
  }

  test_for_validity <- function(index,nc,p) {
    if(p==0) {
      stop("method not applicable to order 0 models!")	  
    }
    if(nc==1) {
      if(sum(index < 1 | index > p) > 0) {
        stop(paste0("index outside the range of lags 1-",p,"!"))   
      }
    } else {
      if(length(index) != 2) {
        stop(paste0("length of index not equal 2!"))	  	
      }
      if(sum(index < 1 | index > nc) > 0) {
        stop(paste0("index outside the range of time series components 1-",nc,"!"))	   
      }
    }
  }
  
  sar_residual <- function(y.qser,fit) {
    p <- fit$p
    n <- fit$n
    ntau <- length(fit$tau)
    if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))  
    nc <- dim(y.qser)[1]  
	if(p > 0) {
      resid <- array(0,dim=c(nc,n,ntau))
      for(l in c(1:ntau)) {
	    if(nc > 1) {
          for(tt in c((p+1):n)) {
            for(j in c(1:p)) resid[,tt,l] <- resid[,tt,l] + fit$A[,,j,l] %*% y.qser[,tt-j,l]
          }
		} else {
          for(tt in c((p+1):n)) {
            for(j in c(1:p)) resid[1,tt,l] <- resid[1,tt,l] + fit$A[j,l] %*% y.qser[1,tt-j,l]
          }		
		}
      }
	  resid <- y.qser - resid
      resid <- resid[,-c(1:p),]
	} else {
	  resid <- y.qser
	}
    resid
  }

  if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))  
  nc <- dim(y.qser)[1]
  n <- dim(y.qser)[2]
  ntau <- dim(y.qser)[3]
  p <- fit$p
  
  test_for_validity(index,nc,p) 
  
  if(method[1]=="ar") {
    # compute residuals from AR model fitted to QSER for each quantile
    resid <- qser2ar(y.qser,p=p)$resid
  } else {
    # compute residuals from SAR model fitted to QSER
    resid <- sar_residual(y.qser,fit)
  }
  if(nc==1) resid <- array(resid,dim=c(1,dim(resid)))
   
  # generate bootstrap sampling of time stamps
  nn <- 2*n
  set.seed(seed)
  idx <- list()
  for(sim.idx in c(1:nsim)) {
    idx[[sim.idx]] <- sample(c(1:dim(resid)[2]),nn,replace=TRUE)
  }
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  if(n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("qser2sar","sar.gc.coef"))
    doParallel::registerDoParallel(cl)
    out <- foreach::foreach(sim.idx=c(1:nsim)) %dopar% {
      if(!mthreads) {
        RhpcBLASctl::blas_set_num_threads(1)
      }
	  qser.sim <- simulate_qser_array(resid,idx[[sim.idx]],fit,index)
      fit.sim <- qser2sar(qser.sim,fit$tau,fit$d,p=fit$p,spar=fit$spar,weighted=fit$weighted)
      A.sim <- sar.gc.coef(fit.sim,index)
    }
    parallel::stopCluster(cl)
  } else {
    out <- foreach::foreach(sim.idx=c(1:nsim)) %do% {
	  qser.sim <- simulate_qser_array(resid,idx[[sim.idx]],fit,index)
      fit.sim <- qser2sar(qser.sim,fit$tau,fit$d,p=fit$p,spar=fit$spar,weighted=fit$weighted)
      A.sim <- sar.gc.coef(fit.sim,index)
    }
  }
  A.sim <- array(0,dim=c(nsim,dim(out[[1]])))
  for(sim.idx in c(1:nsim)) {
    A.sim[sim.idx,,] <- out[[sim.idx]]
  }
  return(A.sim)
}


#' Wald Test and Confidence Band for Granger-Causality Analysis
#'
#' This function computes Wald test and confidence band for Granger-causality 
#' using bootstrap samples generated by \code{sar.gc.bootstrap()} 
#' based the spline autoregression (SAR) model of quantile series (QSER).
#' @param A matrix of selected SAR coefficients 
#' @param A.sim simulated bootstrap samples from \code{sar.gc.bootstrap()}
#' @param sel.lag indices of time lags for Wald test (default = \code{NULL}: all lags)
#' @param sel.tau indices of quantile levels for Wald test (default = \code{NULL}: all quantiles)
#' @return a list with the following elements:
#'   \item{test}{list of Wald test result containing \code{wald} and \code{p.value}}
#'   \item{A.u}{matrix of upper limits of 95\% confidence band of \code{A}}
#'   \item{A.l}{matrix of lower limits of 95\% confidence band of \code{A}}
#' @import 'MASS'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.sar <- qspec.sar(cbind(y1,y2),tau=tau,p=1)
#' A <- sar.gc.coef(y.sar$fit,index=c(1,2))
#' A.sim <- sar.gc.bootstrap(y.sar$qser,y.sar$fit,index=c(1,2),nsim=5)
#' y.gc <- sar.gc.test(A,A.sim)
sar.gc.test <- function(A,A.sim,sel.lag=NULL,sel.tau=NULL) {
   
  wald_test <- function(A,A.sim,sel.lag,sel.tau) {
  # this function computes the bootstrap Wald statistic and p-value
  # of selected elements in A using bootstrap samples in A.sim for Granger Causality test
  # input:   A     = p*ntau-matrix of selected SAR(p) coefficients for all lags and quantiles
  #          A.sim = nsim*p*ntau-array of bootstrap samples associated with A
  #          sel.lap = selected lags for testing (NULL for all lags )
  #          sel.tau = selected quantile levels for testing (NULL for all quantiles)
  # imports: MASS
    if(is.null(sel.lag)) {
      sel.lag <- c(1:dim(A)[1])
    }
    if(is.null(sel.tau)) {
      sel.tau <- c(1:dim(A)[2])
    }
    if(sum(sel.lag < 1 | sel.lag > dim(A)[1]) > 0) {
      stop(paste0("sel.lag outside the range of lag index 1-",dim(A)[1]))	
    }
    if(sum(sel.tau < 1 | sel.tau > dim(A)[2]) > 0) {
      stop(paste0("sel.tau outside the range of quantile index 1-",dim(A)[2]))	
    }
    nlag <- length(sel.lag)
    ntau <- length(sel.tau)
    B <- dim(A.sim)[1]
    Sig <- matrix(0,nrow=ntau*nlag,ncol=ntau*nlag)
    for(b in c(1:B)) {
      Sig <- Sig + tcrossprod(c(A.sim[b,sel.lag,sel.tau]))
    } 
    Sig <- Sig / B
    Sig <- MASS::ginv(Sig)
    chi2 <- c(crossprod(c(A[sel.lag,sel.tau]),Sig) %*% as.vector(c(A[sel.lag,sel.tau])))
    chi2.H0 <- rep(0,B)
    for(b in c(1:B)) {
      chi2.H0[b] <- crossprod(c(A.sim[b,sel.lag,sel.tau]),Sig) %*% as.vector(c(A.sim[b,sel.lag,sel.tau]))
    }
    p.value <- length(which(chi2.H0 >= chi2))/B
    return(list(wald=chi2,p.value=p.value,Sig.inv=Sig,a=c(A[sel.lag,sel.tau])))
  }

  test <- wald_test(A,A.sim,sel.lag,sel.tau)
  A.u <- apply(A.sim,c(2,3),quantile,prob=0.975)
  A.l <- apply(A.sim,c(2,3),quantile,prob=0.025)
  return(list(test=test,A.u=A.u,A.l=A.l))
}


#' Bootstrap Simulation of SAR Coefficients for Testing Equality of Granger-Causality in Two Samples
#'
#' This function simulates bootstrap samples of selected spline autoregression (SAR) coefficients 
#' for testing equality of Granger-causality in two samples based on their SAR models
#' under H0: effect in each sample equals the average effect.
#' @param y.qser matrix or array of QSER from \code{qser()} or \code{qspec.sar()$qser}
#' @param fit object of SAR model from \code{qser2sar()} or \code{qspec.sar()$fit}
#' @param fit2 object of SAR model for the other sample
#' @param index a pair of component indices for multiple time series 
#' or a sequence of lags for single time series (default = \code{c(1,2)})
#' @param nsim number of bootstrap samples (default = 1000)
#' @param method method of residual calculation: \code{"ar"} (default) or \code{"sar"}
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param mthreads if \code{TRUE}, multithread BLAS is enabled when available (default = \code{FALSE}, required for parallel computing) 
#' @param seed seed for random sampling (default = \code{1234567})
#' @return array of simulated bootstrap samples of selected SAR coefficients
#' @import 'foreach'
#' @import 'parallel'
#' @import 'doParallel'
#' @import 'RhpcBLASctl'
#' @export
#' @examples
#' y11 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y21 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y12 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y22 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y1.sar <- qspec.sar(cbind(y11,y21),tau=tau,p=1)
#' y2.sar <- qspec.sar(cbind(y12,y22),tau=tau,p=1)
#' A1.sim <- sar.eq.bootstrap(y1.sar$qser,y1.sar$fit,y2.sar$fit,index=c(1,2),nsim=5)
#' A2.sim <- sar.eq.bootstrap(y2.sar$qser,y2.sar$fit,y1.sar$fit,index=c(1,2),nsim=5)
sar.eq.bootstrap <- function(y.qser,fit,fit2,index=c(1,2),nsim=1000,method=c("ar","sar"),n.cores=1,mthreads=FALSE,seed=1234567) {

  simulate_qser_array_eq <- function(resid,sample.idx,fit,fit2,idx0) {
  # this function simulates quantile series under H0 from bootstrap samples of residuals
  # input: resid = nc*(n-p)*ntau-array of residuals from qser2ar()
  #        sample.idx = nn-vector of sampled time indices
  #        fit, fit2 = SAR models from qser2sar() for sample 1 and sample 2
  #        idx0 = indices of SAR coefficients to set to average under H0
  # output: nc*n*ntau array of simulated quantile series
      p <- fit$p
      A <- fit$A
      A2 <- fit2$A
      nc <- dim(resid)[1]
	  n <- dim(resid)[2] + p
      ntau <- dim(resid)[3]
      nn <- length(sample.idx)
      resid.sim <- array(resid[,sample.idx,],dim=c(nc,nn,ntau))
      qser.sim <- array(0,dim=c(nc,nn,ntau))
      for(ell in c(1:ntau)) {
        for(jj in c((p+1):nn)) {
          if(nc == 1) {
            for(j in c(1:p)) {
              A0 <- A[j,ell]
              if(j %in% idx0) A0 <- 0.5*(A0 + A2[j,ell])
              A0 <- matrix(A0,ncol=nc,nrow=nc)			  
              qser.sim[,jj,ell] <- qser.sim[,jj,ell] + A0 %*% matrix(qser.sim[,jj-j,ell],nrow=nc)
            }		  
          } else {
            for(j in c(1:p)) {
              A0 <- A[,,j,ell]
              A0[idx0[1],idx0[2]] <- 0.5*(A0[idx0[1],idx0[2]] + A2[idx0[1],idx0[2],j,ell])
              qser.sim[,jj,ell] <- qser.sim[,jj,ell] + A0 %*% matrix(qser.sim[,jj-j,ell],nrow=nc)
            } 
          }
          qser.sim[,jj,ell] <- qser.sim[,jj,ell] + resid.sim[,jj,ell] 
        }
      }
      qser.sim <- qser.sim[,c((nn-n+1):nn),]
      return(qser.sim)
  }

  check_for_validity <- function(index,nc,p) {
    if(p==0) {
      stop("method not applicable to order 0 models!")	  
    }
    if(nc==1) {
      if(sum(index < 1 | index > p) > 0) {
        stop(paste0("index outside the range of lags 1-",p,"!"))   
      }
    } else {
      if(length(index) != 2) {
        stop(paste0("length of index not equal 2!"))	  	
      }
      if(sum(index < 1 | index > nc) > 0) {
        stop(paste0("index outside the range of time series components 1-",nc,"!"))	   
      }
    }
  }
    
  sar_residual <- function(y.qser,fit) {
    p <- fit$p
    n <- fit$n
    ntau <- length(fit$tau)
    if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))  
    nc <- dim(y.qser)[1]  
	if(p > 0) {
      resid <- array(0,dim=c(nc,n,ntau))
      for(l in c(1:ntau)) {
	    if(nc > 1) {
          for(tt in c((p+1):n)) {
            for(j in c(1:p)) resid[,tt,l] <- resid[,tt,l] + fit$A[,,j,l] %*% y.qser[,tt-j,l]
          }
		} else {
          for(tt in c((p+1):n)) {
            for(j in c(1:p)) resid[1,tt,l] <- resid[1,tt,l] + fit$A[j,l] %*% y.qser[1,tt-j,l]
          }		
		}
      }
	  resid <- y.qser - resid
      resid <- resid[,-c(1:p),]
	} else {
	  resid <- y.qser
	}
    resid
  }
 
  if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))  
  nc <- dim(y.qser)[1]
  n <- dim(y.qser)[2]
  ntau <- dim(y.qser)[3]
  p <- fit$p
  
  check_for_validity(index,nc,p) 
  
  # compute residauls from sample 1
  if(method[1]=="ar") {
    # from AR model fitted to QSER for each quantile
    resid <- qser2ar(y.qser,p=p)$resid
  } else {
    # from SAR model fitted to QSER
    resid <- sar_residual(y.qser,fit)
  }
  if(nc==1) resid <- array(resid,dim=c(1,dim(resid)))
  
  # generate bootstrap sampling of time stamps
  nn <- 2*n
  set.seed(seed)
  idx <- list()
  for(sim.idx in c(1:nsim)) {
    idx[[sim.idx]] <- sample(c(1:dim(resid)[2]),nn,replace=TRUE)
  }
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  if(n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, c("qser2sar","sar.gc.coef"))
    doParallel::registerDoParallel(cl)
    out <- foreach::foreach(sim.idx=c(1:nsim)) %dopar% {
      if(!mthreads) {
        RhpcBLASctl::blas_set_num_threads(1)
      }
	  qser.sim <- simulate_qser_array_eq(resid,idx[[sim.idx]],fit,fit2,index)
      fit.sim <- qser2sar(qser.sim,fit$tau,fit$d,p=fit$p,spar=fit$spar,weighted=fit$weighted)
      A.sim <- sar.gc.coef(fit.sim,index)
    }
    parallel::stopCluster(cl)
  } else {
    out <- foreach::foreach(sim.idx=c(1:nsim)) %do% {
	  qser.sim <- simulate_qser_array_eq(resid,idx[[sim.idx]],fit,fit2,index)
      fit.sim <- qser2sar(qser.sim,fit$tau,fit$d,p=fit$p,spar=fit$spar,weighted=fit$weighted)
      A.sim <- sar.gc.coef(fit.sim,index)
    }
  }
  A.sim <- array(0,dim=c(nsim,dim(out[[1]])))
  for(sim.idx in c(1:nsim)) {
    A.sim[sim.idx,,] <- out[[sim.idx]]
  }
  return(A.sim)
}

  

#' Wald Test and Confidence Band for Equality of Granger-Causality in Two Samples
#'
#' This function computes Wald test and confidence band for equality of Granger-causality in two samples
#' using bootstrap samples generated by \code{sar.eq.bootstrap()} based on the spline autoregression (SAR) models
#' of quantile series (QSER).
#' @param A1 matrix of selected SAR coefficients for sample 1
#' @param A1.sim simulated bootstrap samples from \code{sar.eq.bootstrap()} for sample 1
#' @param A2 matrix of selected SAR coefficients for sample 2
#' @param A2.sim simulated bootstrap samples from \code{sar.eq.bootstrap()} for sample 2
#' @param sel.lag indices of time lags for Wald test (default = \code{NULL}: all lags)
#' @param sel.tau indices of quantile levels for Wald test (default = \code{NULL}: all quantiles)
#' @return a list with the following elements:
#'   \item{test}{list of Wald test result containing \code{wald} and \code{p.value}}
#'   \item{D.u}{matrix of upper limits of 95\% confidence band for \code{A1 - A2}}
#'   \item{D.l}{matrix of lower limits of 95\% confidence band for \code{A1 - A2}}
#' @import 'MASS'
#' @export
#' @examples
#' y11 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y21 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y12 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y22 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y1.sar <- qspec.sar(cbind(y11,y21),tau=tau,p=1)
#' y2.sar <- qspec.sar(cbind(y12,y22),tau=tau,p=1)
#' A1.sim <- sar.eq.bootstrap(y1.sar$qser,y1.sar$fit,y2.sar$fit,index=c(1,2),nsim=5)
#' A2.sim <- sar.eq.bootstrap(y2.sar$qser,y2.sar$fit,y1.sar$fit,index=c(1,2),nsim=5)
#' A1 <- sar.gc.coef(y1.sar$fit,index=c(1,2))
#' A2 <- sar.gc.coef(y2.sar$fit,index=c(1,2))
#' test <- sar.eq.test(A1,A1.sim,A2,A2.sim,sel.lag=NULL,sel.tau=NULL)
sar.eq.test <- function(A1,A1.sim,A2,A2.sim,sel.lag=NULL,sel.tau=NULL) {
   
  wald_eq_test <- function(A1,A1.sim,A2,A2.sim,sel.lag,sel.tau) {
  # this function computes the two-sample bootstrap Wald statistic and p-value
  # of selected elements in A and A2 using bootstrap samples in A.sim and A2.sim 
  # for testing the equality of A1 and A2, i.e., H0: (A1-A2)[selected elements] = 0
  # input:   A1, A2 = p*ntau-matrix of selected SAR(p) coefficients for all lags and quantiles
  #          A1.sim. A2.sim = nsim*p*ntau-array of bootstrap samples associated with A1 and A2 under H0
  #          sel.lap = selected lags for testing (NULL for all lags )
  #          sel.tau = selected quantile levels for testing (NULL for all quantiles)
  # imports: MASS
    if(is.null(sel.lag)) {
      sel.lag <- c(1:dim(A1)[1])
    }
    if(is.null(sel.tau)) {
      sel.tau <- c(1:dim(A1)[2])
    }
    if(sum(sel.lag < 1 | sel.lag > dim(A1)[1]) > 0) {
      stop(paste0("sel.lag outside the range of lag index 1-",dim(A1)[1]))	
    }
    if(sum(sel.tau < 1 | sel.tau > dim(A1)[2]) > 0) {
      stop(paste0("sel.tau outside the range of quantile index 1-",dim(A1)[2]))	
    }
    nlag <- length(sel.lag)
    ntau <- length(sel.tau)
    B <- dim(A1.sim)[1]
    Sig1 <- matrix(0,nrow=ntau*nlag,ncol=ntau*nlag)
    Sig2 <- matrix(0,nrow=ntau*nlag,ncol=ntau*nlag)
    for(b in c(1:B)) {
      Sig1 <- Sig1 + tcrossprod(c(A1.sim[b,sel.lag,sel.tau]))
      Sig2 <- Sig2 + tcrossprod(c(A2.sim[b,sel.lag,sel.tau]))	  
    } 
    Sig1 <- Sig1 / B
    Sig2 <- Sig2 / B	
    Sig <- MASS::ginv(Sig1+Sig2)
    delta <- c(A1[sel.lag,sel.tau]-A2[sel.lag,sel.tau])
    chi2 <- c(crossprod(delta,Sig) %*% as.vector(delta))
    chi2.H0 <- rep(0,B)
    for(b in c(1:B)) {
	  tmp <- c(A1.sim[b,sel.lag,sel.tau] - A2.sim[b,sel.lag,sel.tau])
      chi2.H0[b] <- crossprod(tmp,Sig) %*% as.vector(tmp)
    }
    p.value <- length(which(chi2.H0 >= chi2))/B
    return(list(wald=chi2,p.value=p.value,Sig.inv=Sig,D=delta))
  }

  test <- wald_eq_test(A1,A1.sim,A2,A2.sim,sel.lag,sel.tau)
  D.u <- apply(A1.sim-A2.sim,c(2,3),quantile,prob=0.975)
  D.l <- apply(A1.sim-A2.sim,c(2,3),quantile,prob=0.025)
  return(list(test=test,D.u=D.u,D.l=D.l))
}


# -- new functions in version 3.0  (October 2024)


#' Periodogram (PER)
#'
#' This function computes the periodogram or periodogram matrix for univariate or multivariate time series.
#' @param y vector (n) or matrix (n x nc) of time series
#' @return A vector (n) or array (nc x nc x n) of periodogram
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y.per <- per(y)
#' plot(y.per)
per <- function(y) {
  if(!is.matrix(y)) y <- matrix(y,ncol=1,nrow=length(y))
  nc <- ncol(y)
  n <- nrow(y)
  y.dft <- matrix(0,nrow=n,ncol=nc)
  for(k in c(1:nc)) {
    y.dft[,k] <- fft(y[,k])
  }  
  y.per <- array(0,dim=c(nc,nc,n))
  for(k in c(1:nc)) {
    for(kk in c(k:nc)) {
	  y.per[k,kk,] <- y.dft[,k] * Conj(y.dft[,kk]) / n
      if(k != kk) y.per[kk,k,] = y.per[k,kk,]
	}
  }
  if(nc == 1) y.per <- Re(y.per[1,1,])
  return(y.per)
}


#' Quantile-Crossing Series (QCSER)
#'
#' This function creates the quantile-crossing series (QCSER) for univariate or multivariate time series.
#' @param y vector or matrix of time series
#' @param tau vector of quantile levels in (0,1)
#' @param normalize \code{TRUE} or \code{FALSE} (default): normalize QCSER to have unit variance
#' @return A matrix or array of quantile-crossing series
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qser <- qcser(y,tau)
#' dim(y.qser)
qcser <- function(y,tau,normalize=FALSE) {
  if(!is.matrix(y)) y <- matrix(y,ncol=1,nrow=length(y))
  n <- nrow(y)
  nc <- ncol(y)
  nq <- length(tau)
  y.qcser <- array(0,dim=c(nc,n,nq))
  for(k in c(1:nc)) {
    q <- quantile(y[,k],tau)
    for(j in c(1:nq)) {
	  y.qcser[k,,j] <- 0
	  y.qcser[k,which(y[,k] <= q[j]),j] <- 1
	  y.qcser[k,,j] <- tau[j] - y.qcser[k,,j]
	  if(normalize) y.qcser[k,,j] <- y.qcser[k,,j] / sqrt(tau[j]*(1-tau[j]))
	}
  }
  if(nc==1) y.qcser <- matrix(y.qcser,ncol=nq)
  return(y.qcser)
}


#' ACF of Quantile Series (QSER) or Quantile-Crossing Series (QCACF)
#'
#' This function creates the ACF of quantile series or quantile-crossing series
#' @param y.qser matrix or array of quantile-crossing series
#' @return A matrix or array of ACF
#' @export
#' @examples
#' y <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' tau <- seq(0.1,0.9,0.05)
#' y.qser <- qcser(y,tau)
#' y.qacf <- qser2qacf(y.qser)
#' dim(y.qacf)
qser2qacf <- function(y.qser) {
  if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,nrow(y.qser),ncol(y.qser)))
  nc <- dim(y.qser)[1]
  n <- dim(y.qser)[2]
  nq <- dim(y.qser)[3]
  y.qacf <- array(0,dim=c(nc,nc,n,nq))
  if(nc > 1) {
    for(j in c(1:nq)) {
      tmp <- stats::acf(t(as.matrix(y.qser[,,j],ncol=nc)),
	                    type="covariance",lag.max=n-1,plot=FALSE,demean=TRUE)$acf
      for(k in c(1:nc)) {
        for(kk in c(1:nc)) {
          y.qacf[k,kk,,j] <- tmp[,k,kk] 
	    }
      }
    }
    rm(tmp)
  } else {
    for(j in c(1:nq)) {
      tmp <- stats::acf(c(y.qser[,,j]),
	                    type="covariance",lag.max=n-1,plot=FALSE,demean=TRUE)$acf 
      y.qacf[1,1,,j] <- tmp[,1,1]
    }	 
  }
  if(nc==1) y.qacf <- matrix(y.qacf[1,1,,],ncol=nq)
  return(y.qacf)
}


# Smoothing of AR Residual Covariance Matrix Across Quantiles
#
# This function perfoms quantile smoothing of Residual Covariance matrix
# @param V array of functional residual covariance to be smoothed
# @param method quantile smoothing method: \code{"gamm"}, \code{"sp"}, or \code{"none"} (default)
Vsmooth <- function(V,method=c("none","gamm","sp")) {
    nc <- dim(V)[1]
	ntau <- dim(V)[3]
	type <- 1
    if(method[1] %in% c("gamm","sp")) {
  	  if(nc == 1) {
	    V[1,1,] <- qsmooth(V[1,1,],method=method[1])
	  } else {
	    if(type[1]==2) {
  	      U <- V
	      for(ell in c(1:ntau)) {
	        tmp <- svd(V[,,ell])
		    U[,,ell] <- tmp$u %*% diag(sqrt(tmp$d))
          }
          for(k in c(1:nc)) {
	        for(kk in c(1:nc)) {   
	          U[k,kk,] <- qsmooth(U[k,kk,],method=method[1])
	        }
	      }
	      for(ell in c(1:ntau)) {
		    V[,,ell] <- U[,,ell] %*% t(U[,,ell])
          }
	    } else {
          for(k in c(1:nc)) {
	        for(kk in c(k:nc)) {   
	          V[k,kk,] <- qsmooth(V[k,kk,],method=method[1])
		      V[kk,k,] <- V[k,kk,]
	        }
	      }	   
	    }
	  }
    }
    return(V)
}	


# Smoothing of AR Coefficient Matrix Across Quantiles
#
# This function perfoms quantile smoothing of AR Coefficients
# @param A array of functional AR coefficients to be smoothed
# @param method quantile smoothing method: \code{"gamm"}, \code{"sp"}, or \code{"none"} (default)
Asmooth <- function(A,method=c("none","gamm","sp")) {
    nc <- dim(A)[1]
    p <- dim(A)[3]
    if(method[1] %in% c("gamm","sp")) {
      for(k in c(1:nc)) {
	    for(kk in c(1:nc)) {
	      for(j in c(1:p)) {
	        A[k,kk,j,] <- qsmooth(A[k,kk,j,],method=method[1])
	      }
        }
	  }
    }
    return(A)
}



#' Autoregression (AR) Model of Quantile Series
#'
#' This function fits an autoregression (AR) model to quantile series (QSER) separately for each quantile level using \code{stats::ar()}.
#' @param y.qser matrix or array of pre-calculated QSER, e.g., using \code{qser()}
#' @param p order of AR model (default = \code{NULL}: selected by AIC)
#' @param order.max maximum order for AIC if \code{p = NULL} (default = \code{NULL}: determined by \code{stats::ar()})
#' @param method quantile smoothing method: \code{"gamm"}, \code{"sp"}, or \code{"NA"} (default)
#' @return a list with the following elements:
#'   \item{A}{matrix or array of AR coefficients}
#'   \item{V}{vector or matrix of residual covariance}
#'   \item{p}{order of AR model}
#'   \item{n}{length of time series}
#'   \item{residuals}{matrix or array of residuals}
#' @import 'stats'
#' @import 'mgcv'
#' @export
qser2ar <- function(y.qser,p=NULL,order.max=NULL,method=c("none","gamm","sp")) {
  if(is.matrix(y.qser)) y.qser <- array(y.qser,dim=c(1,dim(y.qser)))
  nc <- dim(y.qser)[1]
  n <- dim(y.qser)[2]
  ntau <- dim(y.qser)[3]
  # order selection by aic
  if(is.null(p)) {
    aic <- NULL
    for(ell in c(1:ntau)) {
      mu <- apply(matrix(y.qser[,,ell],nrow=nc),1,mean)
      aic <- cbind(aic,stats::ar(t(matrix(y.qser[,,ell]-mu,nrow=nc)),order.max=order.max,method=c("yule-walker"))$aic)
    }
    aic <- apply(aic,1,mean)
    # optimal order minimizes aic
    p <- as.numeric(which.min(aic)) - 1
  }
  # p = 0
  if(p == 0) {
     A <- NULL	
     V <- array(0,dim=c(nc,nc,ntau))
     resid <- array(0,dim=c(nc,n,ntau))
     for(ell in c(1:ntau)) {
       mu <- apply(matrix(y.qser[,,ell],nrow=nc),1,mean)
       V[,,ell] <- matrix(tcrossprod(matrix(y.qser[,,ell]-mu,nrow=nc)),ncol=nc)/n
       resid[,,ell] <- matrix(y.qser[,,ell]-mu,nrow=nc)
     }
     V <- Vsmooth(V,method[1])	 
     if(nc==1) {
       V <- V[1,1,]
       resid <- resid[1,,]
     }
     return(list(A=A,V=V,p=p,n=n,residuals=y.qser))
  }
  # p > 0
  V <- array(0,dim=c(nc,nc,ntau))
  A <- array(0,dim=c(nc,nc,p,ntau))
  resid <- array(0,dim=c(nc,n-p,ntau))
  for(ell in c(1:ntau)) {
    mu <- apply(matrix(y.qser[,,ell],nrow=nc),1,mean)
    fit <- stats::ar(t(matrix(y.qser[,,ell]-mu,nrow=nc)), aic=FALSE, order.max = p, method=c("yule-walker"))
    if(nc==1) {
      V[,,ell] <- sum(fit$resid[-c(1:p)]^2)/n
      for(j in c(1:p)) A[,,j,ell] <- fit$ar[j]	
      resid[,,ell] <- t(fit$resid[-c(1:p)])
    } else {
      V[,,ell] <- matrix(crossprod(fit$resid[-c(1:p),]),ncol=nc)/n
      for(j in c(1:p)) A[,,j,ell] <- fit$ar[j,,]
      resid[,,ell] <- t(fit$resid[-c(1:p),])	  
    }	
  }
  A <- Asmooth(A,method[1])
  V <- Vsmooth(V,method[1])
  if(nc==1) {
    V <- V[1,1,]
    A <- matrix(A[1,1,,],nrow=p)
    resid <- resid[1,,]
  }
  return(list(A=A,V=V,p=p,n=n,residuals=resid))
}


#' Autoregression (AR) Estimator of Quantile Spectrum
#'
#' This function computes autoregression (AR) estimate of quantile spectrum from time series or quantile series (QSER).
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param y.qser matrix or array of pre-calculated QSER (default = \code{NULL}: compute from \code{y} and \code{tau});
#' @param method quantile smoothing method: \code{"gamm"} for \code{mgcv::gamm()}, 
#' \code{"sp"} for \code{stats::smooth.spline()}, or \code{"none"} (default)
#' if \code{y.qser} is supplied, \code{y} and \code{tau} can be left unspecified
#' @param p order of AR model (default = \code{NULL}: automatically selected by AIC)
#' @param order.max maximum order for AIC if \code{p = NULL} (default = \code{NULL}: determined by \code{stats::ar()})
#' @param freq sequence of frequencies in [0,1) (default = \code{NULL}: all Fourier frequencies)
#' @param n.cores number of cores for parallel computing of QDFT if \code{y.qser = NULL} (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing of QDFT (default = \code{NULL})
#' @return a list with the following elements:
#'   \item{spec}{matrix or array of AR quantile spectrum}
#'   \item{freq}{sequence of frequencies}
#'   \item{fit}{object of AR model}
#'   \item{qser}{matrix or array of quantile series if \code{y.qser = NULL}}
#' @import 'stats'
#' @import 'mgcv'
#' @export
#' @examples
#' y1 <- stats::arima.sim(list(order=c(1,0,0), ar=0.5), n=64)
#' y2 <- stats::arima.sim(list(order=c(1,0,0), ar=-0.5), n=64)
#' y <- cbind(y1,y2)
#' tau <- seq(0.1,0.9,0.05)
#' n <- length(y1)
#' ff <- c(0:(n-1))/n
#' sel.f <- which(ff > 0 & ff < 0.5)
#' y.qspec.ar <- qspec.ar(y,tau,p=1)$spec
#' qfa.plot(ff[sel.f],tau,Re(y.qspec.ar[1,1,sel.f,]))
#' y.qser <- qcser(y1,tau)
#' y.qspec.ar <- qspec.ar(y.qser=y.qser,p=1)$spec
#' qfa.plot(ff[sel.f],tau,Re(y.qspec.ar[sel.f,]))
#' y.qspec.arqs <- qspec.ar(y.qser=y.qser,p=1,method="sp")$spec
#' qfa.plot(ff[sel.f],tau,Re(y.qspec.arqs[sel.f,]))
qspec.ar <- function(y,tau,y.qser=NULL,p=NULL,order.max=NULL,freq=NULL,
                method=c("none","gamm","sp"),n.cores=1,cl=NULL) {
  return.qser <- FALSE
  if(is.null(y.qser)) {
    y.qser <- qser(y,tau,n.cores=n.cores,cl=cl)
	return.qser <- TRUE
  }
  fit <- qser2ar(y.qser,p=p,order.max=order.max,method=method[1])
  tmp <- ar2qspec(fit,freq)
  if(return.qser) {
    return(list(spec=tmp$spec,freq=tmp$freq,fit=fit,qser=y.qser)) 
  } else {
    return(list(spec=tmp$spec,freq=tmp$freq,fit=fit))
  }
}



# version 3.1: rename qspec.lw as qacf2speclw
# and absort it into LWQS to become the new qspec.lw  (December 2024)

#' Lag-Window (LW) Estimator of Quantile Spectrum
#'
#' This function computes lag-window (LW) estimate of quantile spectrum
#' with or without quantile smoothing from time series or quantile autocovariance function (QACF).
#' @param y vector or matrix of time series (if matrix, \code{nrow(y)} = length of time series)
#' @param tau sequence of quantile levels in (0,1)
#' @param y.qacf matrix or array of pre-calculated QACF (default = \code{NULL}: compute from \code{y} and \code{tau});
#' if \code{y.qacf} is supplied, \code{y} and \code{tau} can be left unspecified
#' @param M bandwidth parameter of lag window (default = \code{NULL}: quantile periodogram)
#' @param method quantile smoothing method:  \code{"gamm"} for \code{mgcv::gamm()}, 
#' \code{"sp"} for \code{stats::smooth.spline()}, or \code{"none"} (default)
#' @param spar smoothing parameter in \code{smooth.spline()} if \code{method = "sp"} (default = \code{"GCV"})
#' @param n.cores number of cores for parallel computing (default = 1)
#' @param cl pre-existing cluster for repeated parallel computing (default = \code{NULL})
#' @return A list with the following elements:
#'   \item{spec}{matrix or array of spectral estimate}
#'   \item{spec.lw}{matrix or array of spectral estimate without quantile smoothing}
#'   \item{lw}{lag-window sequence}
#'   \item{qacf}{matrix or array of quantile autocovariance function if \code{y.qacf = NULL}}
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
#' y.qacf <- qacf(cbind(y1,y2),tau)
#' y.qper.lw <- qspec.lw(y.qacf=y.qacf,M=5)$spec
#' qfa.plot(ff[sel.f],tau,Re(y.qper.lw[1,1,sel.f,]))
#' y.qper.lwqs <- qspec.lw(y.qacf=y.qacf,M=5,method="sp",spar=0.9)$spec
#' qfa.plot(ff[sel.f],tau,Re(y.qper.lwqs[1,1,sel.f,]))
qspec.lw <- function(y,tau,y.qacf=NULL,M=NULL,method=c("none","gamm","sp"),spar="GCV",n.cores=1,cl=NULL) {
  return.qacf <- FALSE
  if(is.null(y.qacf)) {
    y.qacf <- qacf(y,tau,n.cores=n.cores,cl=cl)
	return.qacf <- TRUE
  }
  y.lw <- qacf2speclw(y.qacf,M)
  if(method[1] %in% c("sp","gamm")) {
    y.qper.lwqs <- qsmooth.qper(y.lw$spec,method=method[1],spar=spar,n.cores=n.cores,cl=cl)
  } else {
    y.qper.lwqs <- y.lw$spec
  }
  if(return.qacf) {
    return(list(spec=y.qper.lwqs,spec.lw=y.lw$spec,lw=y.lw$lw,qacf=y.qacf))
  } else {
    return(list(spec=y.qper.lwqs,spec.lw=y.lw$spec,lw=y.lw$lw))
  }  
}

