#' Composite Likelihood Calculation for Spatial Ordinal Data without Replications (for implmentation)
#'
#' \code{cl_l} Calculate the negative composite log-likelihood value for a particular subject given parameter value and other input variables
#'
#'
#' @param theta a vector of parameter value
#' @param y a vector of observation for the subject.
#' @param X covariate for the particular subject
#' @param dwdv corresponding distance of selected pair
#' @param cmwdv combination of the pairs included into the composite likelihood
#' @param lt number of parameter (i.e. length of theta)
#' @param wn number of pairs
#' @param base identity matrix with dimension \code{J+1}
#' @param J number of category among (ALL) observed response.
#' @param p number of covariate (i.e. number of column of \code{X})
#'
#' @return \code{cl} returns a list: composite log-likelihood value and a vector of first-order partial derivatives for \code{theta}.
#'

cl_l<-function(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p) {

  Xa = X[cmwdv[1,],]
  Xb = X[cmwdv[2,],]
  pairy = matrix(y[cmwdv],nrow=2)
  b = c(-10^4,0,theta[1:(J-2)],10^4,theta[(J-1):(lt-2)])
  Idx=A1lo=A2lo=A1hi=A2hi=list()

  for (j in 0:(J-1)) {
    for (k in 0:(J-1)) {
      A1lo = merge(A1lo,list(base[j+1,]))
      A1hi = merge(A1hi,list(base[j+2,]))
      Idx = merge(Idx, list(which(pairy[1,]==j&pairy[2,]==k)))##among those neighbour pairs
      A2lo = merge(A2lo,list(base[k+1,]))
      A2hi = merge(A2hi, list(base[k+2,]))
    }
  }
  Ldx = rapply(Idx,function(x) length(x),how="list")
  Dwdv2 = rapply(Idx,function(x) dwdv[x],how="unlist")
  X1 = rapply(Idx, function(x) Xa[x,],how="list")
  X2 = rapply(Idx, function(x) Xb[x,],how="list")

  xh <- function(x,y,z) matrix(c(rep(x, each = y),-z), byrow=F,nrow=y)
  X1lo0 = mapply(xh,A1lo,Ldx,X1)
  X1lo <- matrix(NA,0,(J+p+1))
  for( g in 1:(J^2) ){
    if(!length(X1lo0[[g]])==0){
      X1lo <- rbind(X1lo, X1lo0[[g]])
    }
  }
  #X1lo <-do.call(rbind,X1lo0)
  X1hi0 = mapply(xh,A1hi,Ldx,X1)
  X1hi <- matrix(NA,0,(J+p+1))
  for( g in 1:(J^2) ){
    if(!length(X1hi0[[g]])==0){
      X1hi <- rbind(X1hi, X1hi0[[g]])
    }
  }
  #X1hi <- do.call(rbind,X1hi0)
  X2lo0 = mapply(xh,A2lo,Ldx,X2)
  X2lo <- matrix(NA,0,(J+p+1))
  for( g in 1:(J^2) ){
    if(!length(X2lo0[[g]])==0){
      X2lo <- rbind(X2lo, X2lo0[[g]])
    }
  }
  #X2lo <- do.call(rbind,X2lo0)
  X2hi0 = mapply(xh,A2hi,Ldx,X2)
  X2hi <- matrix(NA,0,(J+p+1))
  for( g in 1:(J^2) ){
    if(!length(X2hi0[[g]])==0){
      X2hi <- rbind(X2hi, X2hi0[[g]])
    }
  }
  ### b contains all alfa and beta; alfa0=-10^4; alfaJ = 10^4###
  Mu2hh = cbind(X1hi%*%b,X2hi%*%b) # b contains all alfa and beta
  Mu2ll = cbind(X1lo%*%b,X2lo%*%b) # b contains all alfa and beta
  Mu2hl = cbind(X1hi%*%b,X2lo%*%b)  # b contains all alfa and beta
  Mu2lh = cbind(X1lo%*%b,X2hi%*%b) # b contains all alfa and beta
  rr =  theta[lt]*(1+theta[lt-1]*Dwdv2)*exp(-theta[lt-1]*Dwdv2)
  de = 1-rr^2
  sde = sqrt(de)
  re2 = pbivnorm::pbivnorm(Mu2hh,rho=rr) + pbivnorm::pbivnorm(Mu2ll,rho=rr) -
    (pbivnorm::pbivnorm(Mu2hl,rho=rr)+pbivnorm::pbivnorm(Mu2lh,rho=rr))
  if(sum(is.na(re2))!=0){
    re2[which(is.na(re2))] <- 0
  }
  if(sum(re2)!=0){
    re2[which(re2<=0)]=1e-99 } else re2 <- rep(1e-99, length(re2))  #cl vector value & denominator; length inx

  f2 = -sum(log(re2))/wn

  return(f2) #change
}

