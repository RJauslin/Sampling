### Test for a new implementation of the samplecube method for 
### strongly sparse auxiliary matrices.
### 
### 

rm(list = ls())


#############################################################
#############################################################
#############################################################
#############################################################

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

#############################################################
#############################################################
#############################################################
#############################################################

reduxB <- function(B){
  
  # initialization
  eps <- 1e-12
  sums <- colSums(B)
  sums_row <- rowSums(B)
  B_out <- B
  ind_col <- seq(1,ncol(B_out),1)
  ind_row <- seq(1,nrow(B_out),1)
  
  # loop while any colsums equal to 0 exists
  while(any( sums < eps | sums < -eps)){
    
    # exctract right column
    ind_col <- ind_col[which(sums > eps | sums < -eps)]
    B_out <- B_out[,sums > eps | sums < -eps]
    
    
    # remove duplicated rows
    ind_row <- ind_row[!duplicated(B_out)]
    B_out <- B_out[!duplicated(B_out), ]
    
    # calculate rowsums and remove rowsums equal to 0
    sums_row <- rowSums(B_out)
    ind_row <- ind_row[which(sums_row > eps | sums_row < -eps)]
    
    
    # if we have enough row then compress B, if equal drop one variable then break
    if(length(ind_row) > (ncol(B_out) + 1)){
      B_out <- B_out[which(sums_row > eps | sums_row < -eps)[1:(ncol(B_out) + 1)],]
    }else{
      if(length(ind_row) == length(ind_col)){
        ind_col <- ind_col[-length(ind_col)]
        B_out <- B_out[,-ncol(B_out)]
      }
      break;
    }
    # update sums
    sums <- colSums(B_out)
  }
  # last update rowsums
  sums_row <- rowSums(B_out)
  ind_row <- ind_row[which(sums_row > eps | sums_row < -eps)]
  
  # return the compressed matrix and the column and row that represent B_out.
  return(list(B = B_out,ind_col = ind_col, ind_row = ind_row))
}

#############################################################
#############################################################
#############################################################
#############################################################

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