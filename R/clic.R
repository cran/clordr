#' Composite likelihood Information Criterion
#'
#' \code{clic} Calculating the Composite likelihood information criterion proposed by Varin and Vidoni (2005)
#'
#' @return \code{CLIC}: Composite likelihood information criterion proposed by Varin and Vidoni (2005)
#' @param logCL value of composite log-likelihood.
#' @param mat.hessian hessian matrix.
#' @param mat.J Sensitivity matrix
#'
#' @return \code{clic}: Composite likelihood information criterion proposed by Varin and Vidoni (2005)
#'
#' @details Varin and Vidoni (2005) proposed the information criterion in the form:
#' \eqn{-2*logCL(theta) + 2*trace(H^{-1}(\theta)J(\theta))}
#'
#' @references Varin, C. and Vidoni, P. (2005) A note on composite likelihood inference and model selection. \emph{Biometrika} 92: 519--528.
#'
#' @import pbivnorm  MASS rootSolve parallel doParallel foreach tmvmixnorm utils stats ttutils
#' @export

clic <- function(logCL,mat.hessian,mat.J){
  2*(-logCL + sum(diag(solve(mat.hessian)%*%mat.J)))
}
