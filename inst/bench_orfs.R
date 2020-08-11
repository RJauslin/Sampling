
library(sampling)
library(BalancedSampling)
library(MASS)
library(microbenchmark)

rm(list = ls())
eps=1e-12

# library(devtools)
# install_github("Rjauslin/SamplingC@master")
# library(SamplingC)



#-------------------  INCLUSTION PROBABILITIES

t       <- 3
Pik     <- matrix(rep(0, t*N), ncol = t)
n1 <- N/3
n2 <- N/5
n3 <- N/7
Pik[,1] <- rep(n1/N,N)
Pik[,2] <- inclusionprobabilities(df[,4],n2)
Pik[,3] <- inclusionprobabilities(df[,3],n3)

#-------------------  system.test

system.time(test <- Orfs(Pik,arma = TRUE))
system.time(test <- Orfs(Pik,arma = FALSE))


#-------------------  microbenchmark



