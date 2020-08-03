

rm(list = ls())
library(sampling)
library(MASS)
kk = 1

for(kk in 100:200){
  set.seed(kk)
  eps=0.0000001
  N=1000
  n=300
  p=50
  
  z=runif(N)
  #z=rep(1,N)
  pik=inclusionprobabilities(z,n)
  X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
  A=X/pik
  
  
  test1 <- fastflightcube(X,pik,order=1,comment=TRUE)
  
  if(length(which(test1 > eps & test1 < (1-eps))) < 51){
    print(kk)
    break;
  }
}



set.seed(24)
eps=0.0000001
N=1000
n=300
p=50

z=runif(N)
#z=rep(1,N)
pik=inclusionprobabilities(z,n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik


test1 <- fastflightcube(X,pik,order=1,comment=TRUE)
test2 <- BalancedSampling::flightphase(pik,X)
test3 <- Sampling::flightphase(pik,X)
length(which(test1 > eps & test1 < (1-eps)))
length(which(test2 > eps & test2 < (1-eps)))
length(which(test3 > eps & test3 < (1-eps)))
