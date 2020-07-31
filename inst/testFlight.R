library(BalancedSampling)

#set.seed(12345);
eps = 1e-12
N = 1000; # population size
n = 100; # sample size
pik = rep(n/N,N); # inclusion probabilities
p = 20
X = cbind(pik,matrix(runif(N*p),ncol = p,nrow = N))
pikstar = BalancedSampling::flightphase(pik,X)
pikstar2 = Sampling::flightphase(pik,X)
pikstar3 = BalancedSampling::cube(pik,X)

length(which(pikstar > eps & pikstar < (1-eps)))
length(which(pikstar2 > eps & pikstar2 < (1-eps)))

