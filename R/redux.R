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
#' tmp <- reduxArma(B)
#' tmp2 <- reduxB(B)
#' B_redux <- tmp$B
#' B_redux[1:10,1:20]
#' B[tmp$ind_row[1:10],tmp$ind_col[1:20]]
reduxB <- function(B){
  
  # initialization
  eps <- 1e-12
  sums <- colSums(B)
  sums_row <- rowSums(B)
  B_out <- B
  ind_col <- seq(1,ncol(B_out),1)
  ind_row <- seq(1,nrow(B_out),1)
  
  
  step = 1
  # loop while any colsums equal to 0 exists
  while(any( sums < eps & sums > -eps)){
    # print(step)
    
    # exctract right column
    coltmp <- which(sums > eps | sums < -eps)
    if(length(coltmp)<= 1){
      break;
    }
    B_out <- B_out[,coltmp]
    ind_col <- ind_col[coltmp]

    # calculate rowsums and remove rowsums equal to 0
    sums_row <- rowSums(B_out)
    rowtmp <- which(sums_row > eps | sums_row < -eps)
    B_out <- B_out[rowtmp,]
    ind_row <- ind_row[rowtmp]
    
    # recompute
    sums_row <- rowSums(B_out)
    
    # unique 
    # sums_row <- rowSums(B_out)
    # uniqueRow <- !duplicated(sums_row)
    # ind_row <- ind_row[uniqueRow]
    # B_out <- B_out[uniqueRow, ]
    # sums_row <- rowSums(B_out)
    
    
    
    ## remove duplicated rows
    # uniqueRow <- !duplicated(B_out,MARGIN = 1)
    # ind_row <- ind_row[uniqueRow]
    # B_out <- B_out[uniqueRow, ]
    # sums_row <- rowSums(B_out)
    
    
    
    # if we have enough row then compress B
    if(nrow(B_out) >= (ncol(B_out) + 1)){
      ind_row <- ind_row[1:(ncol(B_out)+1)]
      B_out <- B_out[1:(ncol(B_out)+1),]
    }else{
      if(nrow(B_out) == ncol(B_out)){
        ind_col <- ind_col[-length(ind_col)]
        B_out <- B_out[,-ncol(B_out)]
      }else{
        # ind_row smaller than 
        ind_col <- ind_col[1:(length(ind_row)-1)]
        B_out <- B_out[,1:length(ind_row)-1]
      }
      # break;
    }
    # update sums
    sums <- colSums(B_out)
    step = step + 1;
  }
  # last update rowsums
  # sums_row <- rowSums(B_out)
  # ind_row <- ind_row[which(sums_row > eps | sums_row < -eps)]
  
  # return the compressed matrix and the column and row that represent B_out.
  return(list(B = B_out,ind_col = ind_col, ind_row = ind_row))
}
