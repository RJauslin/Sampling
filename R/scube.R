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
  if (!is.array(X)){
    X = array(X, c(N, length(X)/N))
  }
  p = ncol(X)
  pikstar = pik
  
  
  for (i in 0:(p - 1)) {
    if (length(pikstar[pikstar > EPS & pikstar < (1 - EPS)]) > 0){
      pikstar = ffscube(X[, 1:(p - i)]/pik*pikstar, pikstar, comment)
    }
  }
  pikfin = pikstar
  for (i in 1:N) if (runif(1) < pikfin[i]){
    pikfin[i] = 1
  }
      
  if (comment) {
    A = X[pik > EPS, ]/pik[pik > EPS]
    TOT = t(A) %*% pik[pik > EPS]
    EST = t(A) %*% pikfin[pik > EPS]
    DEV = 100 * (EST - TOT)/TOT
    cat("\n\nQUALITY OF BALANCING\n")
    if (is.null(colnames(X)))
      Vn = as.character(1:length(TOT))
    else Vn = colnames(X)
    for (i in 1:length(TOT)) if (Vn[i] == "")
      Vn[i] = as.character(i)
    d = data.frame(TOTALS = c(TOT), HorvitzThompson_estimators = c(EST),
                   Relative_deviation = c(DEV))
    rownames(d) <- Vn
    print(d)
  }
  round(pikfin,9)
}
