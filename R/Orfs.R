#' @encoding UTF-8
#' @title ORFS method
#'
#' @description
#' Select temporal samples of fixed size at each wave.
#' It allows to generate an optimal temporal rotation in the selection of units using the systematic sampling method.
#'
#'
#' @param Pik a matrix of temporal inclusion probabilities.
#' Columns of \code{Pik} correspond to the waves, and rows correspond to the units.
#' Inclusion probabilities can be totally unequal.
#'
#' @param tol the tolerance parameter. Default value is 1e-9.
#'
#'
#'
#' @return \code{S} a matrix that contains temporal samples.
#' It is the update of \code{Pik} and contains only 0s and 1s that indicates if a unit is selected or not at each wave.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#' @examples
#' \dontrun{
#' ## Temporal inclusion probabilities with 3 waves and 4 units ##
#' Pik <- matrix(c(0.6,0.3,0.3,
#'                 0.2,0.4,0.9,
#'                 0.3,0.2,0.5,
#'                 0.9,0.1,0.3), ncol = 3, byrow = TRUE)
#'
#' N <- 400
#' n1 <- round(N/3)
#' n2 <- round(N/5)
#' n3 <- round(N/7)
#' Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),n1),
#' sampling::inclusionprobabilities(runif(N),n2),
#' sampling::inclusionprobabilities(runif(N),n3)),ncol = 3)
#'
#' ## ORFS method ##
#' system.time(res <- Orfs(Pik, tol = 1e-6))
#' utilisateur     système      écoulé 
#' # 156.53        1.29      157.84
#' # system.time(res <- Orfs(Pik, tol = 1e-6)) #### WITH samplecubeSPOT
#' 
#'}
#' @export
Orfs <- function(Pik, tol = 1e-9)
{
  library(sampling)

  N <- nrow(Pik)
  t <- ncol(Pik)

  # SYSTEMATIC SAMPLING
  res    <- systematicDesign(Pik[1,])
  S      <- as.matrix(res$samples)
  P      <- res$probas
  R      <- rep(1, each = length(res$probas))
  for(i in 2:N){
    res  <- systematicDesign(Pik[i,])
    S    <- rbind(S,res$samples)
    P    <- c(P, res$probas)
    R    <- c(R, rep(i, each = length(res$probas)))
  }
  P      <- as.vector(P)
  R.init <- R
  P.init <- P

  PP     <- TemporalCube(S, P, R)

  return(S[PP == 1,])
}


#' @encoding UTF-8
#' @title ORFS method
#'
#' @export
#' @examples
#' N <- 200
#' Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),20),
#' sampling::inclusionprobabilities(runif(N),30),
#' sampling::inclusionprobabilities(runif(N),40)),ncol = 3)
#' X <- PM(Pik)
PM <- function(Pik, tol = 1e-9)
{
  library(sampling)

  N <- nrow(Pik)
  t <- ncol(Pik)

  # SYSTEMATIC SAMPLING
  res    <- systematicDesign(Pik[1,])
  S      <- as.matrix(res$samples)
  P      <- res$probas
  R      <- rep(1, each = length(res$probas))
  for(i in 2:N){
    res  <- systematicDesign(Pik[i,])
    S    <- rbind(S,res$samples)
    P    <- c(P, res$probas)
    R    <- c(R, rep(i, each = length(res$probas)))
  }
  P      <- as.vector(P)

  RR <- unique(R)
  Z  <- matrix(rep(0,nrow(S)*length(RR)), nrow = nrow(S))

  for(i in 1:length(RR)){
    Z[R == RR[i],i] <- 1
  }

  M           <- cbind(Z,S)
  colnames(M) <- NULL


  return(list(P = P,PM = P*M))
}
