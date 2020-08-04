
#' Title
#'
#' @param X
#' @param pik
#' @param order
#' @param comment
#'
#' @return
#' @export
#'
#' @examples
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
fastflightcubeSPOT <- function (X, pik, order = 1, comment = TRUE)
{
  EPS = 1e-11

  "reduc" <- function(X) {
    EPS = 1e-10
    N = dim(X)[1]
    Re = svd(X)
    array(Re$u[, (Re$d > EPS)], c(N, sum(as.integer(Re$d >
                                                      EPS))))
  }


  ########################### START ALGO


  N = length(pik)
  p = round(length(X)/length(pik))
  X <- array(X, c(N, p))
  if (order == 1){
    o <- sample(N, N)
  }else {
    if (order == 2){
      o <- seq(1, N, 1)
    }else {
      o <- order(pik, decreasing = TRUE)
    }
  }
  liste <- o[(pik[o] > EPS & pik[o] < (1 - EPS))]
  if (comment == TRUE) {
    cat("\nBEGINNING OF THE FLIGHT PHASE\n")
    cat("The matrix of balanced variable has", p, " variables and ",
        N, " units\n")
    cat("The size of the inclusion probability vector is ",
        length(pik), "\n")
    cat("The sum of the inclusion probability vector is ",
        sum(pik), "\n")
    cat("The inclusion probability vector has ", length(liste),
        " non-integer elements\n")
  }
  pikbon <- pik[liste]
  Nbon = length(pikbon)
  Xbon <- array(X[liste, ], c(Nbon, p))
  pikstar <- pik
  flag = 0


  # X <- Xbon
  # pik <- pikbon



  # begin algorithm (general case where N > p) at the end of this phase you have
  # at most p values that are not equal to 0 or 1.
  if (Nbon > p) {
    if (comment == TRUE)
      cat("Step 1  ")
    system.time(pikstarbon <- algofastflightcube2(Xbon, pikbon))
    # system.time(pikstarbon <- algofastflightcube(Xbon, pikbon))
    pikstar[liste] = pikstarbon
    flag = 1
  }

  # reupdate the liste and the exctract element no equal to 0 or 1
  liste <- o[(pikstar[o] > EPS & pikstar[o] < (1 - EPS))]
  pikbon <- pikstar[liste]
  Nbon = length(pikbon)
  Xbon <- array(X[liste, ], c(Nbon, p))
  pbon = dim(Xbon)[2]

  # if you still have value that are not equal to 0 or 1 you reduc the matrix and loop until
  if (Nbon > 0) {
    Xbon = reduc(Xbon)
    pbon = dim(Xbon)[2]
  }
  k = 2
  while (Nbon > pbon & Nbon > 0) {
    if (comment == TRUE)
      cat("Step ", k, ",  ")
    k = k + 1
    pikstarbon <- algofastflightcube2(Xbon/pik[liste] * pikbon,
                                     pikbon)
    pikstar[liste] = pikstarbon
    liste <- o[(pikstar[o] > EPS & pikstar[o] < (1 - EPS))]
    pikbon <- pikstar[liste]
    Nbon = length(pikbon)
    Xbon <- array(X[liste, ], c(Nbon, p))
    if (Nbon > 0) {
      Xbon = reduc(Xbon)
      pbon = dim(Xbon)[2]
    }
    flag = 1
  }
  if (comment == TRUE)
    if (flag == 0)
      cat("NO FLIGHT PHASE")
  if (comment == TRUE)
    cat("\n")
  pikstar
}
