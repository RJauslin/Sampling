#' Title
#'
#' @param X
#' @param pik
#' @param deepness
#' @param EPS
#'
#' @return
#' @export
#'
#' @examples
#'
#' rm(list = ls())
#' library(sampling)
#' library(MASS)
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
#'
#' Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
#' B=cbind(Z,-Z)
#' r=c(ceiling(pik%*%B))
#' r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
#' s=fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.000001)
#' round(s,3)
#' c(t(B)%*%pik)
#' c(r)
#' c(t(B)%*%s)
#' colSums(X)
#' c(t(X)%*%(s/pik))
#'
#' piks=fast.flight.cube(cbind(X,B*pik),pik,deepness=1,EPS=0.000001)
#' pikstar=fast.flight.cube.ineq(X/pik*piks,piks,B,r,deepness=1,EPS=0.000001)
fast.flight.cube<-function(X,pik,deepness=1,EPS=0.0000001)
{
  A=as.matrix(X/pik)
  if(nrow(A)==0){
    A=matrix(0,c(length(pik),1))
  }
  if(is.null(X)){
    prof=1
  }else{
    if(is.matrix(X)){
      prof=ncol(X)+1
    }else{
        prof=2
    }
  }
  TEST=(EPS<pik) & (pik<1-EPS)
  prof2=min(sum(TEST),prof)
  if(prof2==0){
    a=0
  }else{
    pikR=pik[TEST][1:prof2]
    AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,]
    NN=MASS::Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))))
    a=ncol(NN)
  }

  while(a>0.5){
    u=NN[,1] # extract column of kernel

    l1=min(pmax((1-pikR)/u,-pikR/u))
    l2=min(pmax((pikR-1)/u,pikR/u))
    if(runif(1)<l2/(l1+l2)){
      pik[TEST][1:prof2] = pikR+l1*u
    }else{
      pik[TEST][1:prof2] = pikR-l2*u
    }


    TEST=(EPS<pik) & (pik<1-EPS) # redefine 0 entries
    prof2=min(sum(TEST),prof) # redefine dimension of submatrix B

    if(prof2==0){
      a=0
    }else{
      pikR=pik[TEST][1:prof2]
      AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,]
      NN=MASS::Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))))
      a=ncol(NN)
      # print(a)
    }
  }
  pik
}
