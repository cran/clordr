#' Surrogate Residuals for Replciations of Spatial Ordinal Data
#'
#' \code{surrogate.residual} simulate the surrogate residual with the the given parameter value and covariate for model diagnostics.
#'
#' @param response a matrix of observation (row: spatial site and column: subject).
#' @param covar regression (design) matrix, including intercepts.
#' @param location a matrix contains spatial location of sites within each subject.
#' @param seed seed input for simulation (default =\code{NULL}).
#' Parameter values:
#' @param midalpha cutoff for latent ordinal response.
#' @param beta regression coefficient  for \code{covar}.
#'
#' @param sigma2 \eqn{\sigma^2} for exponential covariance.
#' @param phi spatial correlation for exponential covariance.
#' @param burn.in burn-in length (i.e. declaring the initial sample).
#' @param output logical flag indicates whether printing out result (default: \code{TRUE}).
#'
#' @details Given vector of observed responses, the design matrix, spatial location for sites and parameter value, raw surrogate residuals are simulated using an efficient Gibbs sampling, which can be used for model diagnostics. When the fitted model is correct, the raw surrogate residuals among subjects should follow multivariate normal with mean 0 and covariance Sigma. If the model is correct, residual plot should be close to a null plot or random scatter. For example, it can be used to check the potential missing in covariate, non-linearity of covariate and outliers. In particular for the example below, the residual plot shows that linearity of Xi is adequate for the model.
#'
#' @return \code{surrogate.residual} returns a (no. spatial site * no. subject) matrix contains
#' raw surrogate residuals with element corresponds to the response matrix.
#'
#' @import pbivnorm  MASS rootSolve parallel doParallel foreach tmvmixnorm utils stats ttutils
#' @export
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
#' response <- sim.rord(n.subject, n.site, n.rep = 1,
#' midalpha, beta, phi, sigma2, covar=VV, location)[[1]]
#'
#' \donttest{
#' # Example for linearity of covariate
#' sur.resid <- surrogate.residual(response, covar=VV, location, seed =1,
#' midalpha, beta, sigma2, phi,
#' burn.in=20, output = TRUE)
#'
#' scatter.smooth(rep(Xi,each=n.site),c(sur.resid),
#' main="Surrogate residual against Xi", xlab="Xi", ylab="Surrogate residual",
#' lpars = list(col = "red", lwd = 3, lty = 2))
#'
#' abline(h=0, col="blue")
#' }
#'

surrogate.residual <- function(response, covar, location, seed = NULL, midalpha, beta, sigma2, phi,
                               burn.in=20, output = TRUE) {

  # n: n.subject
  # N: number of sites for each subject
  n <- ncol(response) ; N <- nrow(response)

  residual <- matrix(NA, nrow = N, ncol = n)

  talpha1 <- c(-Inf, 0, midalpha, Inf)
  J <- length(talpha1)
  start.value <-  (talpha1[1:(J-1)] + talpha1[2:J])/2 ;  start.value[1] <- talpha1[2] - 0.1 ; start.value[J-1] <-  talpha1[J-1] + 0.1


  L2matrix = as.matrix(dist(location, method='euclidean', upper=TRUE, diag=TRUE))
  SIGMA = sigma2*(1+phi*L2matrix)*exp(-phi*L2matrix)

  set.seed(seed)

  for(i in 1:n){
    y.tmp <- response[,i] ; X.tmp <- covar[((i-1)*N+1):(i*N),]

    lower.tmp <- talpha1[y.tmp +1] ; upper.tmp <- talpha1[y.tmp + 2]
    start.tmp <- start.value[y.tmp+1]



    residual[,i] <- tmvmixnorm::rtmvn(n=1, Mean=c(X.tmp%*%beta), Sigma=SIGMA, D = diag(1, N),
                                      lower=lower.tmp, upper=upper.tmp, int=start.tmp,burn = burn.in, thin = 0) - c(X.tmp%*%beta)

      #c(tmvtnorm::rtmvnorm(n=1, mean=c(X.tmp%*%beta), sigma=SIGMA,
      #                                   lower=lower.tmp, upper=upper.tmp, algorithm="gibbs",
      #                                   start.value=start.tmp, burn.in.samples = burn.in)) - c(X.tmp%*%beta)

    if (output) {
      cat("Surrogate residual generated:",i,"out of", n,"\r")
      flush.console()
    }

  }

  if (output) cat("\b\n")

  residual
}


