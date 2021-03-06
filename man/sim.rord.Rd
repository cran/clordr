% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.rord.R
\name{sim.rord}
\alias{sim.rord}
\title{Simulation of Replciations of Spatial Ordinal Data}
\usage{
sim.rord(
  n.subject,
  n.site,
  n.rep = 100,
  midalpha,
  beta,
  phi,
  sigma2,
  covar,
  location
)
}
\arguments{
\item{n.subject}{number of subjects.}

\item{n.site}{number of spatial sites for each subject.}

\item{n.rep}{number of simulation.
Parameter inputs include:}

\item{midalpha}{cutoff parameter (excluding -Inf and +Inf);}

\item{beta}{regression coefficient;}

\item{phi}{dependence parameter for spatial dependence; and}

\item{sigma2}{sigma^2 (variance) for the spatial dependence.}

\item{covar}{regression (design) matrix, including intercepts.}

\item{location}{a matrix contains spatial location of sites within each subject.}
}
\value{
\code{sim.rord} returns a list (length \code{n.rep}) of matrix (\code{n.subject*n.site}) with the underlying parameter as inputs.
}
\description{
\code{sim.rord} Simulate replications of spatial ordinal data
}
\examples{
set.seed(1203)
n.subject <- 100
n.lat <- n.lon <- 10
n.site <- n.lat*n.lon

beta <- c(1,2,-1) # First 1 here is the intercept
midalpha <- c(1.15, 2.18) ; phi <- 0.8 ; sigma2 <- 0.7

true <- c(midalpha,beta,sigma2,phi)

Xi <- rnorm(n.subject,0,1) ; Xj <- rbinom(n.site,1,0.6)

 VV <- matrix(NA, nrow = n.subject*n.site, ncol = 3)

 for(i in 1:n.subject){ for(j in 1:n.site){
     VV[(i-1)*n.site+j,] <- c(1,Xi[i],Xj[j])
       }
 }

location <- cbind(rep(seq(1,n.lat,length=n.lat),n.lat),rep(1:n.lon, each=n.lon))
sim.data <- sim.rord(n.subject, n.site, n.rep = 2, midalpha, beta, phi, sigma2, covar=VV, location)

length(sim.data)
head(sim.data[[1]])
dim(sim.data[[1]])
hist(sim.data[[1]])
}
