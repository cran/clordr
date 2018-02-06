#' Composite Likelihood Estimation for Replciations of Spatial Ordinal Data
#'
#' \code{cle.rord} Estimate parameters (including regression coefficient and cutoff) for replications of spatial ordinal data using pairwise likelihood approach. Initial estimate is calculated by assuming ther is no spatial dependence (using \code{MASS::polr()})
#'
#' @param response a matrix of observation (row: spatial site and column: subject).
#' @param covar regression (design) matrix, including intercepts.
#' @param location a matrix contains spatial location of sites within each subject
#' @param radius radius for selecting pairs for the composite likelihood estimation.
#' @param n.sim number of simulation used for parametric bootstrapping (and hence used for asymptotic variance and standard error).
#' @param output logical flag indicates whether printing out result (default: \code{TRUE}).
#'
#' @details Given vector of ordinal responses, the design matrix, spatial location for sites, weight radius (for pair selection), and the prespecified number of simulation used for estimating the Godambe information matrix. The function first estimates parameters of interest by maximizing the composite log-likelihood using \code{optim(...,method = "L-BFGS-B")}, then computes the simulated based standard error and asymptotic covariance matrix by parametric boostrapping.
#' @return \code{cle.rord} returns a list contains:
#' @return \code{vec.par}: a vector of estimator for \eqn{\theta=(\alpha,\beta,\sigma^2,\phi)};
#' @return \code{vec.se}: a vector of standard error for the estimator;
#' @return \code{mat.asyvar}: estimated asymptotic covariance matrix \eqn{H^{-1}(\theta)J(\theta)H^{-1}(\theta)} for the estimator; and
#' @return \code{mat.Hessian}: Hessian matrix at the parameter estimate.
#' @return \code{mat.J}: Sensitivity matrix estimated by parametric boostrapping.
#' @return \code{CLIC}: Composite likelihood information criterion (see help manual of \code{clic()} for detail)
#'
#' @examples
#' set.seed(1203)
#' n.subject <- 50
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
#' \dontrun{
#' options(digits=3)
#' result <- cle.rord(response=sim.data[[1]], covar=VV,
#'           location ,radius = 4, n.sim = 100, output = FALSE)
#' result$vec.par
#' # alpha2  alpha3   beta0   beta1   beta2 sigma^2     phi
#' # 1.195   2.223   0.880   1.859  -0.956   0.754   0.828
#'
#' result$vec.se
#' # alpha2  alpha3   beta0   beta1   beta2 sigma^2     phi
#' # 0.0469  0.0677  0.0700  0.0884  0.0446  0.1024  0.0557
#'
#' }
#'
#'
cle.rord <- function(response, covar, location ,radius=4, n.sim=100, output = TRUE) {

  n <- ncol(response) ; N <- nrow(response) ; J <- length(levels(factor(response))) ; p <- NCOL(covar)
  lt <- J + p

  i <- j <- 1

  d <- c(dist(location, method='euclidean')) # calculate the distance between each pair of sites

  wd <- which(d <= radius) # select pairs of sites to be included
  dwdv <- t(d[wd]) # corresponding distance of selected pair
  wn <- length(dwdv) # number of pairs
  cm <- combn(1:N,2)   # the combination of the pairs
  cmwdv <- cm[,wd] # the combination of the pairs included into the composite likelihood
  base <- diag(J+1)   # the identity matrix

  # Initial value
  po1 = MASS::polr(factor(as.vector(response) )~covar[,-1],method='probit')
  ini <- unname(c(po1$zeta[-1] - po1$zeta[1],-po1$zeta[1],po1$coef,0.762,0.826))


  func<-function(theta, response, covar){
    sum_func <- 0
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_func <- sum_func+cl_l(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)
    }
    return(sum_func)
  }


  dfun<-function(theta, response, covar){
    sum_gr <- numeric(J+p)
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_gr <- sum_gr+cl(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)$gr
    }
    return(sum_gr)
  }


  opt <- optim(par=ini, fn=func, gr=dfun, lower=ini-abs(ini)*0.5,
        upper=ini+abs(ini)*0.1,
        method = 'L-BFGS-B', response=response, covar=covar)

  vec.par <- unname(opt$par)

  if (output) {
    cat("Estimate for\n")
    cat("alpha:",round(vec.par[1:(J-2)],digits=3),"\n")
    cat("beta:",round(vec.par[(J-1):(J+p-2)],digits=3),"\n")
    cat("sigma^2, phi:",round(tail(vec.par,2),digits=3),"\n")

  }

  sim.tmp <- sim.rord(n.subject=n, n.site=N, n.rep=n.sim,
                      midalpha=vec.par[1:(J-2)], beta=vec.par[(J-1):(lt-2)], sigma2=vec.par[lt-1], phi=vec.par[lt],
                      covar, location)

  score.tmp <- matrix(0, ncol=lt, nrow=n.sim)

  for (j in 1:n.sim) {
    response.tmp <- sim.tmp[[j]]
    score.tmp[j,] <- numDeriv::grad(function(x){func(theta=x,response=response.tmp, covar)}, x=vec.par,method="simple")
  }

  mat.H <- numDeriv::hessian(function(x){func(theta=x,response=response, covar)}, vec.par)
  mat.H.inv <- solve(mat.H)

  mat.J <- cov(score.tmp)

  mat.asyvar <- mat.H.inv %*% mat.J %*% mat.H.inv

  vec.se <- sqrt(diag(mat.asyvar))

  names(vec.se) <- names(vec.par) <- colnames(mat.asyvar) <- c(
    paste0(rep("alpha",J-2),2:(J-1)),paste0(rep("beta",p),0:(p-1)),"sigma^2","phi")

  CLIC <- clic(logCL= -opt$value,mat.hessian=mat.H,mat.J=mat.J)


list(vec.par=vec.par, vec.se=vec.se, mat.asyvar=mat.asyvar, mat.Hessian=mat.H, mat.J=mat.J, CLIC=CLIC)
}
