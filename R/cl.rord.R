#' Composite Likelihood Calculation for Replciations of Spatial Ordinal Data (for illustration)
#'
#' \code{cl.rord} Calculate the  composite log-likelihood value for replications of spatial ordinal data at given value of parameter value.
#' Note that this function is not directly used in \code{cle.rord} but illustration only.
#'
#' @param theta a vector of parameter value
#' @param response a matrix of observation (row: spatial site and column: subject).
#' @param covar regression (design) matrix, including intercepts.
#' @param location a matrix contains spatial location of sites within each subject
#' @param radius radius for selecting pairs for the composite likelihood estimation.
#'
#' @return \code{cl.rord} returns a list: negative composite log-likelihood, a vector of first-order partial derivatives for \code{theta}.
#'
#' @examples
#' set.seed(1203)
#' n.subject <- 10
#' n.lat <- n.lon <- 10
#' n.site <- n.lat*n.lon
#'
#' beta <- c(1,2,-1) # First 1 here is the intercept
#' midalpha <- c(1.15, 2.18) ; sigma2 <- 0.7 ; phi <- 0.8
#'
#' true = c(midalpha,beta,sigma2,phi)
#'
#' Xi = rnorm(n.subject,0,1) ; Xj <- rbinom(n.site,1,0.6)
#'
#'  VV <- matrix(NA, nrow = n.subject*n.site, ncol = 3)
#'
#'  for(i in 1:n.subject){ for(j in 1:n.site){
#'      VV[(i-1)*n.site+j,] <- c(1,Xi[i],Xj[j])
#'        }
#'  }
#'
#' location = cbind(rep(seq(1,n.lat,length=n.lat),n.lat),rep(1:n.lon, each=n.lon))
#' sim.data <- sim.rord(n.subject, n.site, n.rep = 2, midalpha, beta, sigma2, phi, covar=VV, location)
#'
#' cl.rord(theta=true,response=sim.data[[1]], covar=VV, location, radius = 4)
#'

cl.rord <- function(theta,response, covar, location ,radius = 4) {

  n <- ncol(response) ; N <- nrow(response) ; J <- length(levels(factor(response))) ; p <- NCOL(covar)

  lt <- J + p

  d <- c(dist(location, method='euclidean')) # calculate the distance between each pair of sites

  wd <- which(d <= radius) # select pairs of sites to be included
  dwdv <- t(d[wd]) # corresponding distance of selected pair
  wn <- length(dwdv) # number of pairs
  cm <- combn(1:N,2)   # the combination of the pairs
  cmwdv <- cm[,wd] # the combination of the pairs included into the composite likelihood
  base <- diag(J+1)   # the identity matrix

  # the negative composite log-likelihood value
  func<-function(theta, response, covar){
    sum_func <- 0
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_func <- sum_func+cl_l(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)
    }
    return(sum_func)
  }

  # the derivative of the negative composite log-likelihood value
  dfun<-function(theta, response, covar){
    sum_gr <- numeric(J+p)
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_gr <- sum_gr+cl(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)$gr
    }
    return(unname(sum_gr))
  }

  sum_func <- func(theta, response, covar)
  sum_gr <- dfun(theta, response, covar)

list(func=sum_func,score=sum_gr)
}

