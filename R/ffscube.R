#' Fast flight phase for the cube method modified
#'
#' @description 
#' 
#' implementation modified from the package sampling.
#' 
#' @param X matrix of auxiliary variables.
#' @param pik vector of inclusion probabilities.
#' @param order order to rearrange the data. Default 1
#' @param comment bool, if comment should be written.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' # Matrix of balancing variables
#' X=cbind(c(1,1,1,1,1,1,1,1,1),c(1,2,3,4,5,6,7,8,9))
#' # Vector of inclusion probabilities.
#' # The sample size is 3.
#' pik=c(1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3)
#' # pikstar is almost a balanced sample and is ready for the landing phase
#' pikstar=fastflightcube(X,pik,order=1,comment=TRUE)
#' round(pikstar,9)
#'
#'
#'
#'
#'
#' rm(list = ls())
#' set.seed(1)
#' eps <- 1e-13
#' library(Matrix)
#' N <- 300
#' Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
#' sampling::inclusionprobabilities(runif(N),50),
#' sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
#' X <- PM(Pik)$PM
#' pik <- PM(Pik)$P
#' dim(X)
#' order = 2
#' EPS = 1e-11
#' 
#' system.time(test1 <- fastflightcubeSPOT(X,pik,order = 2))
#' system.time(test2 <- fastflightcubeSPOT(X,pik,order = 1))
#' system.time(test3 <- sampling::fastflightcube(X,pik, order = 2))
#' system.time(test4 <- BalancedSampling::flightphase(pik,X))
#'
#'
#' }
ffscube <- function (X, pik, comment = TRUE)
{
  
  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  
  X <- as.matrix(X)
  EPS = 1e-11
  N = length(pik)
  p = ncol(X)
  X <- array(X, c(N, p))
  #---- order
  o <- seq(1, N, 1)
  
  
  liste <- o[(pik[o] > EPS & pik[o] < (1 - EPS))]
  pikbon <- pik[liste]
  Nbon = length(pikbon)
  Xbon <- array(X[liste, ], c(Nbon, p))
  pikstar <- pik
  
  
  
  ##---------------------------------------------------------------
  ##                    Redux for the N-p value                   -
  ##---------------------------------------------------------------
  
  
  if(Nbon > p){
    pikstarbon <- algofastflightcubeSPOT(Xbon, pikbon,redux = TRUE)
    pikstar[liste] = pikstarbon
  }
  
  
  ##----------------------------------------------------------------
  ##                            update                             -
  ##----------------------------------------------------------------
  
  
  liste <- o[(pikstar[o] > EPS & pikstar[o] < (1 - EPS))]
  pikbon <- pikstar[liste]
  Nbon = length(pikbon)
  Xbon <- array(X[liste, ], c(Nbon, p))
  pbon = dim(Xbon)[2]
  
  
  ##----------------------------------------------------------------
  ##                Finish flightphase without redux               -
  ##----------------------------------------------------------------
  
  
  while (Nbon > pbon & Nbon > 0) {
    pikstarbon <- algofastflightcubeSPOT(Xbon/pik[liste] * pikbon,pikbon,redux = FALSE)
    pikstar[liste] = pikstarbon
    liste <- o[(pikstar[o] > EPS & pikstar[o] < (1 - EPS))]
    pikbon <- pikstar[liste]
    Nbon = length(pikbon)
    Xbon <- array(X[liste, ], c(Nbon, p))
    pbon = ncol(Xbon)
  }
  
  
  return(pikstar)
}
