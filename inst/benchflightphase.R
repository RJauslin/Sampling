
library(sampling)
library(BalancedSampling)
library(MASS)
library(microbenchmark)
library(Matrix)

rm(list = ls())
eps=1e-12

library(devtools)
remove.packages("SamplingC")
install_github("Rjauslin/SamplingC@master")
library(SamplingC)


#-------------------  
#-------------------  set up 1
#-------------------  
rm(list = ls())
N = 500
n = 80
p = 100
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))



#-------------------  system.time

system.time(test1 <- BalancedSampling::flightphase(pik,X))
system.time(test2 <- SamplingC::ffphase(pik,X))
system.time(test3 <- SamplingC::flightphase_arma(X,pik))
system.time(test4 <- fastflightcube(X,pik))



set.seed(12345);
N = 100000; # population size
n = 100; # sample size
pik = rep(n/N,N); # inclusion probabilities
# matrix of 5 auxiliary variables
X = cbind(p,runif(N),runif(N),runif(N),runif(N)); 

system.time(test1 <- BalancedSampling::flightphase(pik,X))
system.time(test2 <- SamplingC::ffphase(pik,X))
system.time(test4 <- fastflightcube(X,pik))


#-------------------  microbenchmark

microbenchmark(
  BalancedSampling::flightphase(pik,X),
  ffphase(pik,X),
  flightphase_arma(X,pik)
  # fastflightcubeSPOT(X,pik,order = 2,comment = FALSE)
  # fastflight cube(X,pik)
)


#-------------------  
#-------------------  set up 1
#-------------------  


rm(list = ls())
N <- 300
n1 <- round(N/3,N)
n2 <- round(N/5,N)
n3 <- round(N/7,N)

Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),n1),
sampling::inclusionprobabilities(runif(N),n2),
sampling::inclusionprobabilities(runif(N),n3)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P

#-------------------  system.time

system.time(test1 <- BalancedSampling::flightphase(pik,X))
system.time(test2 <- SamplingC::ffphase(pik,X,order = TRUE,redux = FALSE))
system.time(test3 <- SamplingC::flightphase_arma(X,pik))
system.time(test4 <- fastflightcubeSPOT(X,pik,order = 2,comment = FALSE))
system.time(test5 <- sampling::fastflightcube(X,pik,order = 2,comment = FALSE))

A <- X/pik
t(A)%*%pik
t(A)%*%test4

microbenchmark(
  BalancedSampling::flightphase(pik,X),
  ffphase(pik,X),
  flightphase_arma(X,pik),
  fastflightcubeSPOT(X,pik,order = 2,comment = FALSE,method = 2)
  # sampling::fastflightcube(X,pik)
)

