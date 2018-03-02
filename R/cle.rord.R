#' Composite Likelihood Estimation for Replciations of Spatial Ordinal Data
#'
#' \code{cle.rord} Estimate parameters (including regression coefficient and cutoff) for replications of spatial ordinal data using pairwise likelihood approach.
#'
#' @param response a matrix of observation (row: spatial site and column: subject).
#' @param covar regression (design) matrix, including intercepts.
#' @param ini.sp initial estimate for spatial parameter, \eqn{\phi,\sigma^2} (default: c(0.5,0.5)).
#' @param location a matrix contains spatial location of sites within each subject.
#' @param radius radius for selecting pairs for the composite likelihood estimation.
#' @param n.sim number of simulation used for parametric bootstrapping (and hence used for asymptotic variance and standard error).
#' @param output logical flag indicates whether printing out result (default: \code{TRUE}).
#' @param SE logical flag for detailed output.
#' @param parallel logical flag indicates using parallel processing (default: \code{TRUE}).
#' @param n.core number of physical cores used for parallel processing (when \code{parallel} is \code{TRUE}, default value is \code{max(detectCores()/2,1)}).
#' @param est.method logical flag (default) \code{TRUE} for rootsolve and \code{FALSE} for L-BFGS-B.
#' @param maxiter maximum number of iterations in the root solving of gradient function (dafault: 100).
#' @param rtol relative error tolerrance  in the root solving of gradient function (default: 1e-6).
#' @param factr reduction in the objective (-logCL) within this factor of the machine tolerance for L-BFGS-B (default: 1e7).

#'
#' @details Given vector of ordinal responses, the design matrix, spatial location for sites, weight radius (for pair selection), and the prespecified number of simulation used for estimating the Godambe information matrix. Initial estimate is obtained by fitting model without spatial dependence (using \code{MASS::polr()}) and optional guess of spatial parameters. The function first estimates parameters of interest by either solving the gradient of composite log-likelihood using \code{rootSolve::multiroot()} or maximize the composite log-likelihood by \code{optim(..., method="L-BFGS-B")}. The asymptotic covariance matrix and standard error of parameters are then estimated by parametric boostrapping. Although the default root solving option is typically more efficient, it may encounter runtime error if negative value of \eqn{\phi} is evaluated (and L-BFGS-B approach should be used).
#' @return \code{cle.rord} returns a list contains:
#' @return \code{vec.par}: a vector of estimator for \eqn{\theta=(\alpha,\beta,\phi,\sigma^2)};
#'
#' @return \code{vec.se}: a vector of standard error for the estimator;
#' @return \code{mat.asyvar}: estimated asymptotic covariance matrix \eqn{H^{-1}(\theta)J(\theta)H^{-1}(\theta)} for the estimator;
#' @return \code{mat.Hessian}: Hessian matrix at the parameter estimate;
#' @return \code{mat.J}: Sensitivity matrix estimated by parametric boostrapping; and
#' @return \code{CLIC}: Composite likelihood information criterion (see help manual of \code{clic()} for detail).

#'
#' @examples
#' set.seed(1228)
#' n.subject <- 50
#' n.lat <- n.lon <- 10
#' n.site <- n.lat*n.lon
#'
#' beta <- c(1,2,-1) # First 1 here is the intercept
#' midalpha <- c(1.15, 2.18) ; phi <- 0.6 ; sigma2 <- 0.7
#'
#' true <- c(midalpha,beta,phi,sigma2)
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
#' sim.data <- sim.rord(n.subject, n.site, n.rep = 2, midalpha, beta, phi, sigma2, covar=VV, location)
#'
#' \dontrun{
#' options(digits=3)
#' result <- cle.rord(response=sim.data[[1]], covar=VV,
#'           location ,radius = 4, n.sim = 100, output = TRUE, n.core = 4)
#' result$vec.par
#' # alpha2  alpha3   beta0   beta1   beta2     phi sigma^2
#' # 1.175   2.217   1.028   2.036  -1.053   0.586   0.681
#'
#' result$vec.se
#' # alpha2  alpha3   beta0   beta1   beta2     phi sigma^2
#' # 0.0514  0.0873  0.0963  0.1144  0.0462  0.0281  0.0620
#'
#' }
#'

cle.rord <- function(response, covar, location ,radius=4, n.sim=100, output = TRUE, SE = TRUE,
                     parallel = TRUE, n.core = max(detectCores()/2,1),
                     ini.sp = c(0.5,0.5), est.method = TRUE, maxiter = 100, rtol = 1e-6, factr = 1e7) {

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
  ini <- unname(c(po1$zeta[-1] - po1$zeta[1],-po1$zeta[1],po1$coef,ini.sp))


  func<-function(theta, response, covar){

    sum_func <- 0
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_func <- sum_func+cl_l(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)
    }
    return(sum_func)
  }


  dfun <- function(theta, response, covar){

    if(theta[lt-1] < 0) stop("phi cannot be negative, please use est.method = FALSE for L-BFGS-B and/or change ini.sp!!!\n")

    sum_gr <- numeric(J+p)
    for(i in 1:n){
      y <- response[,i]
      X <- covar[((i-1)*N+1):(i*N),]

      sum_gr <- sum_gr+cl(theta,y,X,dwdv,cmwdv,lt,wn,base,J,p)$gr
    }

    return(sum_gr)
  }

  if (est.method) {
  vec.par <- unname(rootSolve::multiroot(f=dfun,start=ini,response=response,covar=covar,
                                         maxiter = maxiter, rtol = rtol)$root)
  } else {
    lower <- ini-abs(ini)*0.5
    upper <- ini+abs(ini)*0.2

    if (lower[lt-1]  < 0) lower[lt-1] <- 1e-4

    vec.par <- unname( optim(par=ini, fn=func, gr=dfun, lower=lower, upper=upper,
                 method = 'L-BFGS-B',
                 control = list(factr = factr),
                 response=response, covar=covar)$par )
  }


  if (output) {
    cat("Estimate for\n")
    cat("alpha:",round(vec.par[1:(J-2)],digits=3),"\n")
    cat("beta:",round(vec.par[(J-1):(J+p-2)],digits=3),"\n")
    cat("phi, sigma^2:",round(tail(vec.par,2),digits=3),"\n")

  }
  if (SE) {
  sim.tmp <- sim.rord(n.subject=n, n.site=N, n.rep=n.sim,
                      midalpha=vec.par[1:(J-2)], beta=vec.par[(J-1):(lt-2)], phi=vec.par[lt-1], sigma2=vec.par[lt],
                      covar, location)

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {(abs(x - round(x)) < tol)&&(!(x== 0))}
  # Small function to check if positive integer

  if (parallel == TRUE){
    i <- 1
    options(warn=-1) # Remove warning message: closing unused connection

    if (!is.double(n.core)) {
      cat("Wrong input type for n.core, replaced by the default value, max(detectCores()/2,1) = ",max(parallel::detectCores()/2,1),"\n")
      n.core <- max(parallel::detectCores()/2,1)
    }else if(!is.wholenumber(n.core)) {
      cat("Input for n.core is not a positive integer, replaced by the default value, max(detectCores()/2,1) = ",max(parallel::detectCores()/2,1),"\n")
      n.core <- max(parallel::detectCores()/2,1)
    }

    cl <- parallel::makeCluster(n.core)
    doParallel::registerDoParallel(cl)
    mat.J <- cov(foreach(i = 1:n.sim, .export="cl", .combine=rbind) %dopar% {
      dfun(theta=vec.par, response=sim.tmp[[i]],covar=covar)
    })

    stopCluster(cl)

    options(warn=0) # Resume warning message

  } else{

  score.tmp <- matrix(0, ncol=lt, nrow=n.sim)

  for (j in 1:n.sim) {
    response.tmp <- sim.tmp[[j]]
    score.tmp[j,] <- dfun(theta = vec.par, response = response.tmp, covar = covar)
  }

  mat.J <- cov(score.tmp)

  }

  if (output) cat("Completed calculation of variability matrix, Now start to calculate Hessian matrix.\n")

  func.outsq <- function(vec) vec %o% vec

  mat.H <- matrix(0,nrow=length(vec.par), ncol=length(vec.par))

  for (i in 1:n) {
    y <- response[,i]
    X <- covar[((i-1)*N+1):(i*N),]

    mat.score.i <- cl(theta=vec.par,y,X,dwdv,cmwdv,lt,wn,base,J,p)$grmat[,c(3:J,(J+2):(lt+3))]


    mat.H <- mat.H + matrix(rowSums(apply(mat.score.i, 1, func.outsq)), nrow = length(vec.par),
                                        byrow = TRUE)/wn

  }


  mat.H.inv <- solve(mat.H)



  mat.asyvar <- mat.H.inv %*% mat.J %*% mat.H.inv

  vec.se <- sqrt(diag(mat.asyvar))
  CLIC <- clic(logCL= func(theta=vec.par,response = response, covar=covar),
               mat.hessian = mat.H,mat.J = mat.J)
}

  names(vec.par) <- c(paste0(rep("alpha",J-2),2:(J-1)),paste0(rep("beta",p),0:(p-1)),"phi", "sigma^2")

if (SE) {
  names(vec.se) <- colnames(mat.asyvar) <- names(vec.par)

    ans <- list(vec.par=vec.par, vec.se=vec.se, mat.asyvar=mat.asyvar, mat.Hessian=mat.H, mat.J=mat.J, CLIC=CLIC)
  } else {

    ans <-  list(vec.par=vec.par)
  }
ans
}
