
library(sampling)
library(BalancedSampling)
library(MASS)
library(microbenchmark)

rm(list = ls())
eps=1e-12


library(devtools)
install_github("Rjauslin/SamplingC@master")
library(SamplingC)


##############
#
# Some examples to show that the function flightphase of package BalancedSampling
# does not goes until the kernel is empty. We propose an alternative function that
# have approximately the same computational time.
# 

set.seed(1)
N <-  1000
n <-  300
p <-  1
q <-  7
z <-  runif(N)

pik <-  inclusionprobabilities(z,n)
X <-  cbind(pik,matrix(rnorm(N*p),c(N,p)))
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
X <- cbind(X,B*pik)


piks <- round(fastflightcube(X,pik),9)
piksBal <- round(BalancedSampling::flightphase(pik,X),9)
piksSam <- round(SamplingC::flightphase(pik,X),9)
pikCube <- round(BalancedSampling::cube(pik,X),9)


dim(X)
length(which(piks > eps & piks < (1-eps)))
length(which(piksBal > eps & piksBal < (1-eps)))
length(which(piksSam > eps & piksSam < (1-eps)))





##############
#
# Microbenchmark



set.seed(1)
N <-  3000
n <-  300
p <-  30
q <-  14
z <-  runif(N)

pik <-  inclusionprobabilities(z,n)
X <-  cbind(pik,matrix(rnorm(N*p),c(N,p)))
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
X <- cbind(X,B*pik)





system.time(piks <- round(fastflightcube(X,pik),9))
system.time(piksBal <- round(BalancedSampling::flightphase(pik,X),9))
system.time(piksSam <- round(SamplingC::flightphase(pik,X),9))
system.time(test <- flightphase_arma(X,pik))
# system.time(pikCube <- round(BalancedSampling::cube(pik,X),9))


dim(X)
length(which(piks > eps & piks < (1-eps)))
length(which(piksBal > eps & piksBal < (1-eps)))
length(which(piksSam > eps & piksSam < (1-eps)))
length(which(test > eps & test < (1-eps)))




