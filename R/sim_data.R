#' Simulation of Replciations of Spatial Ordinal Data
#'
#' \code{sim.rord} Simulate replications of spatial ordinal data
#'
#' @param n.subject number of subjects.
#' @param n.site number of spatial sites for each subject.
#' @param n.rep number of simulation.
#' Parameter inputs include:
#' @param midalpha cutoff parameter (excluding -Inf and +Inf);
#' @param beta regression coefficient;
#' @param sigma2 sigma^2 (variance) for the spatial dependence; and
#' @param phi dependence parameter for spatial dependence.
#' @param covar regression (design) matrix, including intercepts.
#' @param location a matrix contains spatial location of sites within each subject.
#'
#' @return \code{sim.rord} returns a list (length \code{n.rep}) of matrix (\code{n.subject*n.site}) with the underlying parameter as inputs.
#'
#' @examples
#' set.seed(1203)
#' n.subject <- 100
#' n.lat <- n.lon <- 10
#' n.site <- n.lat*n.lon
#'
#' beta <- c(1,2,-1) # First 1 here is the intercept
#' midalpha <- c(1.15, 2.18) ; sigma2 <- 0.7 ; phi <- 0.8
#'
#' true <- c(midalpha,beta,sigma2,phi)
#'
#' Xi <- rnorm(n.subject,0,1) ; Xj <- rbinom(n.site,1,0.6)
#'
#'  VV <- matrix(NA, nrow = n.subject*n.site, ncol = 3)
#'
#'  for(i in 1:n.subject){ for(j in 1:n.site){
#'      VV[(i-1)*n.site+j,] <- c(1,Xi[i],Xj[j])
#'        }
#'  }
#'
#' location <- cbind(rep(seq(1,n.lat,length=n.lat),n.lat),rep(1:n.lon, each=n.lon))
#' sim.data <- sim.rord(n.subject, n.site, n.rep = 2, midalpha, beta, sigma2, phi, covar=VV, location)
#'
#' length(sim.data)
#' head(sim.data[[1]])
#' dim(sim.data[[1]])
#' hist(sim.data[[1]])



sim.rord <- function(n.subject, n.site, n.rep = 100, midalpha, beta, sigma2, phi, covar, location) {

  N <- n.site ; n <- n.subject

  observe<-function(x,alpha) {
    b = rep(0,length(x))
    I = length(alpha)
    alpha=c(0,alpha)
    for(i in 1:I){
      b[which(x>=alpha[i] & x<alpha[i+1])]=i
    }
    b[which(x>=alpha[i+1])]=I+1
    return(b)
  }

  talpha1 = c(-Inf, 0, midalpha, Inf)
  true = c(midalpha,beta,phi,sigma2)
  lt = length(true)


  L2matrix = as.matrix(dist(location, method='euclidean', upper=TRUE, diag=TRUE))
  SIGMA = sigma2*(1+phi*L2matrix)*exp(-phi*L2matrix)


  p <- NCOL(covar) ; J <- length(talpha1)

  R <- chol(SIGMA)

  ystar_gen <- function(x, R=R, beta=beta, sigma2=sigma2){
    zstar <-  t(R)%*%array(rnorm(N),c(N,1))+ covar[((x-1)*N+1):(x*N),]%*%beta
    epstar <- diag(sqrt(1-sigma2),N)%*%array(rnorm(N),c(N,1))
    return(ystar = epstar+zstar)
  }

  result <- NULL

  i <- j <- 1

  for (j in 1:n.rep) {

  ystar <- matrix(0, N, n)

    for(i in 1:n){ ystar[,i] <- ystar_gen(i, R=R, beta=beta, sigma2=sigma2)  }

    yn <- apply(ystar,2,observe,midalpha) # Convert to ordinal

    result <- c(result,list(yn))
  }

  result
}
