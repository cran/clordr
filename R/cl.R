#' Composite Likelihood Calculation for Spatial Ordinal Data without Replications (for implmentation)
#'
#' \code{cl} Calculate the  composite log-likelihood value and score function for a particular subject given parameter value and other input variables.
#'
#' @param theta a vector of parameter value.
#' @param y a vector of observation for the subject.
#' @param X covariate for the particular subject.
#' @param dwdv corresponding distance of selected pair.
#' @param cmwdv combination of the pairs included into the composite likelihood.
#' @param lt number of parameter (i.e. length of theta).
#' @param wn number of pairs with distance.
#' @param base identity matrix with dimension \code{J+1}.
#' @param J number of category among (ALL) observed response.
#' @param p number of covariate (i.e. number of column of \code{X}).
#'
#' @return \code{cl} returns a list: composite log-likelihood value and a vector of first-order partial derivatives for \code{theta}.
#'


cl<-function(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p) {

  xi<-function(u1,u2,r,sqr) (u2-r*u1)/sqr

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
  Mu2hl= cbind(X1hi%*%b,X2lo%*%b)  # b contains all alfa and beta
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

  numer12 = (c(dnorm(Mu2hh[,1])*pnorm(xi(Mu2hh[,1],Mu2hh[,2], rr, sde)))*X1hi+
               c(dnorm(Mu2hh[,2])*pnorm(xi(Mu2hh[,2],Mu2hh[,1], rr, sde)))*X2hi) +
    (c(dnorm(Mu2ll[,1])*pnorm(xi(Mu2ll[,1],Mu2ll[,2], rr, sde)))*X1lo+
       c(dnorm(Mu2ll[,2])*pnorm(xi(Mu2ll[,2],Mu2ll[,1], rr, sde)))*X2lo) -
    (c(dnorm(Mu2hl[,1])*pnorm(xi(Mu2hl[,1],Mu2hl[,2], rr, sde)))*X1hi+
       c(dnorm(Mu2hl[,2])*pnorm(xi(Mu2hl[,2],Mu2hl[,1], rr, sde)))*X2lo) -
    (c(dnorm(Mu2lh[,1])*pnorm(xi(Mu2lh[,1],Mu2lh[,2], rr, sde)))*X1lo+
       c(dnorm(Mu2lh[,2])*pnorm(xi(Mu2lh[,2],Mu2lh[,1], rr, sde)))*X2hi)

  dense = (1/(2*pi*sde))*{
    exp(-0.5*(rowSums(Mu2hh^2)-2*rr*Mu2hh[,1]*Mu2hh[,2])/de) +
      exp(-0.5*(rowSums(Mu2ll^2)-2*rr*Mu2ll[,1]*Mu2ll[,2])/de) -
      exp(-0.5*(rowSums(Mu2hl^2)-2*rr*Mu2hl[,1]*Mu2hl[,2])/de) -
      exp(-0.5*(rowSums(Mu2lh^2)-2*rr*Mu2lh[,1]*Mu2lh[,2])/de)
  }
  numer3 =  c(dense*(1+theta[lt-1]*Dwdv2)*exp(-theta[lt-1]*Dwdv2) )
  numer4 = c(-dense*theta[lt]*Dwdv2*Dwdv2*theta[lt-1]*exp(-theta[lt-1]*Dwdv2))
  grmat2 = cbind(numer12,numer3,numer4)/re2
  f2 = -sum(log(re2))/wn
  gr2 = -colSums(grmat2[,c(3:J,(J+2):(lt+3))])/wn
  return(list(func = f2, gr = gr2, grmat = grmat2)) #change
}

