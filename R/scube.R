#' @title Sample cube SPOT
#'
#' @param X matrix of auxiliary variable
#' @param pik vector of inclusion probabilities
#' @param order order 
#' @param comment comment if true
#' @param method method
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' N=1000
#' p=7
#' set.seed(1)
#' x=rnorm(N*p,10,3)
#' # random inclusion probabilities
#' pik= runif(N)
#' X=array(x,c(N,p))
#' X=cbind(cbind(X,rep(1,times=N)),pik)
#' pikfin = sampling::samplecube(X,pik,1,TRUE)
#' pikfin2 = samplecubeSPOT(X,pik,1,TRUE)
#'
#'
#' rm(list = ls())
#' set.seed(1)
#' eps <- 1e-13
#' library(Matrix)
#' N <- 200
#' Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
#' sampling::inclusionprobabilities(runif(N),50),
#' sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
#' X <- PM(Pik)$PM
#' pik <- PM(Pik)$P
#' dim(X)
#' order = 2
#' method = 1
#' comment = TRUE
#' EPS = 1e-11
#'
#' system.time(s <- samplecubeSPOT(X,pik,order = 2,method = 2,arma = FALSE))
#' system.time(s <- sampling::samplecube(X,pik,order = 2,comment = TRUE,method = 2))
#' system.time(s <- BalancedSampling::cube(pik,X))
#' }
scube <- function(X, pik,comment = FALSE)
{
  
  ##----------------------------------------------------------------
  ##                        initialization                         -
  ##----------------------------------------------------------------
  
  
  EPS = 1e-11
  N = length(pik)
  X = as.matrix(X)
  p = ncol(X)
  pikstar = pik
  
  
  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------
  
  for (i in 0:(p - 1)) {
    if (length(pikstar[pikstar > EPS & pikstar < (1 - EPS)]) > 0){
      pikstar = ffscube(X[, 1:(p - i)]/pik*pikstar, pikstar, comment)
    }
  }
  
  
  ##---------------------------------------------------------------
  ##                      Put equal to 0 or 1                     -
  ##---------------------------------------------------------------
  
  
  pikstar[pikstar < EPS] <- 0
  pikstar[pikstar > (1-EPS)] <- 1
  
  
  return(pikstar)
}
