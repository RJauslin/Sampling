
library(sampling)
library(BalancedSampling)
library(MASS)
library(microbenchmark)

rm(list = ls())
eps=1e-12


library(devtools)
install_github("Rjauslin/SamplingC@master")
library(SamplingC)



rm(list = ls())
N = 500
n = 80
p = 100
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
system.time(test1 <- BalancedSampling::flightphase(pik,X))
system.time(test2 <- flightphaseArma(pik,X))
system.time(test2 <- flightphase_arma(X,pik))
system.time(test2 <- fastflightcube(X,pik))
# 


microbenchmark(
  BalancedSampling::flightphase(pik,X),
  flightphaseArma(pik,X),
  flightphase_arma(X,pik)
  # fastflight cube(X,pik)
)


A <- X/pik

t(A)%*%pik
t(A)%*%test1
t(A)%*%test2
t(A)%*%pik # correc
