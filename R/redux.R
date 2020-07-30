#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
#' rm(list = ls())
#' set.seed(1)
#' eps <- 1e-13
#' library(Matrix)
#' N <- 50
#' Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),5),
#' sampling::inclusionprobabilities(runif(N),10),
#' sampling::inclusionprobabilities(runif(N),15)),ncol = 3)
#' X <- PM(Pik)$PM
#' pik <- PM(Pik)$P
#' dim(X)
#' order = 2
#' EPS = 1e-11
#'
#'
#' p <- ncol(X)
#' A <-  X/pik
#' B <- A[1:(p + 1), ]
#'
#' tmp <- reduxB(B)
#' B_redux <- tmp$B
#' B_redux[1:10,1:20]
#' B[tmp$ind_row[1:10],tmp$ind_col[1:20]]
reduxB <- function(B){
  eps <- 1e-12
  sums <- colSums(B)
  sums_row <- rowSums(B)
  B_out <- B
  ind_col <- seq(1,ncol(B_out),1)
  ind_row <- seq(1,nrow(B_out),1)
  while(any( sums < eps | sums < -eps)){

    ind_col <- ind_col[which(sums > eps | sums < -eps)]
    B_out <- B_out[,sums > eps | sums < -eps]


    # remove duplicated rows
    ind_row <- ind_row[!duplicated(B_out)]
    B_out <- B_out[!duplicated(B_out), ]
    sums_row <- rowSums(B_out)
    ind_row <- ind_row[which(sums_row > eps | sums_row < -eps)]

    if(length(ind_row) > (ncol(B_out) + 1)){
      B_out <- B_out[which(sums_row > eps | sums_row < -eps)[1:(ncol(B_out) + 1)],]
    }else{
      if(length(ind_row) == length(ind_col)){
        ind_col <- ind_col[-length(ind_col)]
        B_out <- B_out[,-ncol(B_out)]
      }
      # print("break")
      break;
    }
    # print(dim(B_out))
    sums <- colSums(B_out)
  }
  sums_row <- rowSums(B_out)
  ind_row <- ind_row[which(sums_row > eps | sums_row < -eps)]
  return(list(B = B_out,ind_col = ind_col, ind_row = ind_row))
}
