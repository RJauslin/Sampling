#' @title Cube method for temporal samples selection
#'
#' @description
#' Select one temporal sample for each unit using the Cube method.
#' The function \code{\link[sampling:samplecube]{samplecube}} from the package \code{sampling} is used.
#'
#'
#' @param S a matrix that contains variables on which the sample must be balanced. See Details.
#'
#' @param P a vector of probabilities of select samples.
#'
#' @param R a vector that specify to which unit each sample belongs.
#'
#' @param tol the tolerance parameter. Default value is 1e-8.
#'
#'
#' @details
#' Balancing constraints considered in \code{S} allow first to select exactly one temporal sample for each unit considered
#' and then to verify constraint of fixed samples size.
#'
#'
#' @return \code{PP} the balanced sample selected (a vector of 0s and 1s with the same size as \code{P}).
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#' @examples
#' \dontrun{
#' 
#' ## Temporal samples of three units ##
#' S <- matrix(c(0,1,1,
#'               1,1,0,
#'               1,0,1,
#'               0,0,1,
#'               0,1,0,
#'               1,0,1,), ncol=3, byrow = T)
#' R <- c(1,1,2,2,3,3)
#' P <- c(0.2,0.8,0.6,0.4,0.5,0.5)
#'
#' ## Find balanced sample ##
#' res  <- LandTemporalPivot(S, P, R, tol = 1e-6)
#' res
#' }
#'
#' @export

TemporalCube <- function(S, P, R, tol = 1e-8){
  library(sampling)
  
  RR <- unique(R)
  Z  <- matrix(rep(0,nrow(S)*length(RR)), nrow = nrow(S))
  
  for(i in 1:length(RR)){
    Z[R == RR[i],i] <- 1
  }
  
  M           <- cbind(Z,S)
  colnames(M) <- NULL
  
  PP               <- samplecubeSPOT(P*M, P,order = 2, method = 2, comment = F)
  # PP               <- sampling::samplecube(P*M, P, order = 2,method = 2, comment = F)
  PP[PP < tol]     <- 0
  PP[PP > (1-tol)] <- 1
  return(PP)
}
