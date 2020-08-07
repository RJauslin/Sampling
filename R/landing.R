#' Lading Ineq (concept function for now)
#'
#' @param X matrix of auxiliary variables
#' @param pikstar vector of inclusion probabilities
#' @param pik vector of inclusion probabilities
#' @param B matrix
#' @param r vector
#' @param comment bool, should comment be written.
#'
#' @return sample (vector of 0 or 1)
#' @export
#'
#' @examples
#' \dontrun{
#' rm(list = ls())
#' EPS=0.0000001
#' N=1000
#' n=300
#' p=1
#' q=7
#' z=runif(N)
#' #z=rep(1,N)
#' pik=inclusionprobabilities(z,n)
#' X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
#' A=X/pik
#' Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
#' B=cbind(Z,-Z)
#' r=c(ceiling(pik%*%B))
#' r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
#' # s=as.vector(ineq(X,pik,B,r,EPS=0.000001))
#'
#' t(B)%*%pik <= r
#' t(B)%*%s <= r
#' colSums(X)
#' c(t(X)%*%(s/pik))
#'
#'
#' # ON COMMENCE PAR EQUILIBRER TOUTEs
#' # LES VARIABLES (Y COMPRIS LES INEGALITE)
#' # A la fin de la phase de vol on a donc au maxium p+q+1 unités non égal 0 ou 1
#'
#' piks <- as.vector(flightphase_arma2(cbind(X,B*pik),pik))
#' p+q+1 # + 1 pour la colonne 1
#' length(which(piks > EPS & piks <(1-EPS)))
#' # check
#' t(B)%*%piks <= r
#' t(X/pik)%*%piks
#' t(X/pik)%*%pik
#'
#' # ON RELACHE LES CONTRAINTES D'INEGALITE
#' # il reste au maxium p+1 unités non égal à 0 ou 1
#' pikstar <- ineq(X/pik*piks,piks,B,r)
#' p+1 # + 1 pour la colonne 1
#' length(which(pikstar > EPS & pikstar <(1-EPS)))
#' # check
#' t(B)%*%pikstar <= r
#' t(X/pik)%*%pikstar
#' t(X/pik)%*%pik
#'
#' s <- landIneq(X,pikstar,pik,B,r)
#' t(B)%*%s <= r
#' t(X/pik)%*%s
#' t(X/pik)%*%pik
#'}
landIneq <- function(X, pikstar, pik,B,r,comment = TRUE){
  # Idea of the algorithm is to first find all sample that have size of sum(pikstar)
  # on element of pikstar that does not have value equal 1 or 0. -> pikland
  #
  # Then we list all the sample that have this size with the function samplen
  #
  # Finally we check the sample that satisfy the inequality i.e. t(B)%*%s <= r.
  #
  # -> firstly study the behaviour
  # if no sample are found TODO
  #
  
  
  # initializing and error
  
  EPS = 1e-11
  p = dim(X)[2]
  N = dim(X)[1]
  q = dim(B)[2]
  if(N != dim(B)[1]){
    stop("The size of the matrix B and X are not the same.")
  }
  index = (pikstar > EPS & pikstar < (1 - EPS))
  index01 =  which(index != TRUE)
  
  pikland = pikstar[index]
  Nland = length(pikland)
  if(comment){
    cat("Number of non-equal 0 or 1 units ", Nland, "\n")
  }
  
  Xland = array(X[index, ], c(Nland, p))
  Bland = array(B[index, ], c(Nland, q))
  rland = r - t(B[index01,])%*%pikstar[index01]
  nland = sum(pikland)
  
  # calculate sampleSet and sampleSetSize
  FLAGI = (abs(nland - round(nland)) < EPS)
  if (FLAGI) {
    pikland = round(nland) * pikland/nland
    nland = round(nland)
    # SSS = writesample(nland, Nland)
    sampleSet = samplen(Nland,nland)
  }else{
    sampleSet = cbind(samplen(Nland,trunc(nland)), samplen(Nland,trunc(nland) + 1))
  }
  sampleSetSize = ncol(sampleSet)
  
  if (comment) {
    cat("The linear program can considered ", sampleSetSize, " possible samples\n")
  }
  
  sampleSettmp <- c()
  for(i in 1:sampleSetSize){
    if(all(t(Bland)%*%sampleSet[,i] <= rland)){
      sampleSettmp <- cbind(sampleSettmp,sampleSet[,i])
    }
  }
  if(is.null(sampleSettmp)){
    # stop("No sample that satisfy constraint")
    message("No sample that satisfy constraint \n Constrained are relaxed")
    print(sampleSetSize)
  }else{
    sampleSet <- sampleSettmp
    sampleSetSize = ncol(sampleSet)
  }
  if (comment) {
    cat("The linear program will considered ", sampleSetSize, " possible samples that satisfies the inequalities\n")
  }
  
  
  # calculate cost
  Astar = matrix(0, p, sampleSetSize)
  Astar <-  t(Xland/pik[index]) %*%(sampleSet-pikland)
  A = X[index,]/pik[index]
  # A = X[pik > EPS, ]/pik[pik > EPS]
  cost = apply(Astar,
               MARGIN = 2,
               FUN = function(x,H){return(t(x)%*%H%*%x)},
               H =  MASS::ginv(t(A) %*% A))
  
  # cost =
  
  # find solution to the system ->  add column of 1 to the fixed sampe size ?
  # find minimum x such that t(V)%"%x <= b
  
  # V = rbind(sampleSet, rep(1, times = sampleSetSize))
  # b = c(pikland, 1)
  # constdir = rep("==", times = (Nland + 1))
  V = sampleSet
  b = pikland
  # solve(V,b)
  
  # Vsvd <- svd(V)
  # Vdiag <- diag(1/Vsvd$d)
  # x <- Vsvd$v %*% Vdiag %*% t(Vsvd$u) %*% b
  # V%*%x
  # b
  
  #
  # V = t(rbind(sampleSet, rep(1, times = sampleSetSize)))
  # b = c(pikland, 1)
  constdir = rep("==", times = (Nland))
  # x = lpSolve::lp("min", rep(1,length(cost)), V, constdir, b)
  x = lpSolve::lp("min", cost, V, constdir, b)
  if(x$status == 2){
    stop("Error: no feasible solution reached.")
  }else{
    x <- x$solution
  }
  if (comment) {
    cat(V%*%x,"\n")
    cat(b,"\n")
  }
  
  # Vtmp <- xtmp <- c()
  # for(i in 1:ncol(V)){
  #   if(all(t(Bland)%*%V[,i] <= rland)){
  #     Vtmp <- cbind(Vtmp,V[,i])
  #     xtmp <- c(xtmp,x[i])
  #   }
  # }
  # if (comment) {
  #   cat(xtmp,"\n")
  # }
  
  # choose a sampleSet randomly
  u = runif(1, 0, 1)
  # cumsum(x)
  # cumsum(x)> u
  i = 0
  ccc = 0
  while (ccc < u) {
    i = i + 1
    ccc = ccc + x[i]
  }
  
  
  # update value
  pikfin = pikstar
  pikfin[index] = sampleSet[,i]
  return(pikfin)
  
  
  
  t(B)%*%pikfin <= r
  t(X/pik)%*%pikfin
  t(X/pik)%*%pik
}

