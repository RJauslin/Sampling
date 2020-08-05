### Test for a new implementation of the samplecube method for 
### strongly sparse auxiliary matrices.
### 
### 

rm(list = ls())
library(Rcpp)

#############################################################
#############################################################
#############################################################
#############################################################

cppFunction(depends = "RcppArmadillo",code = 'Rcpp::List systematicDesign(arma::vec pik){

  //Value to round value like -1e-17...
  double eps = 1e-13;
  double R = 1e+9;

  
  
  
  //index of pik equal 1 or 0
  arma::uvec ones = arma::find(pik > (1-eps));
  arma::uvec zeros = arma::find(pik < eps);

  // add difference between ceiling(pik) and sum(pik) (possibly 0)
  arma::vec pik1(1);
  pik1.fill(ceil(arma::sum(pik)) - sum(pik));
  arma::vec piks = join_cols(pik, pik1);

  // Number of element
  int N = piks.n_elem;

  // vk and r
  arma::vec vk = arma::cumsum(piks);


  //round value otherwise -1e-17 can appears
  vk = round(R*vk)/R;
  arma::vec vk1 = vk - arma::floor(vk);
  arma::vec r = arma::sort(vk1);

  // add 1 to r
  arma::vec one(1);
  one.fill(1);
  r = join_cols(r, one);

  // centered value and p (rounded p to be sure that is correctly equal to 0)
  arma::vec cent = (r(arma::span(0, N-1)) + r(arma::span(1, N)))/2;
  arma::vec p = r(arma::span(1, N)) - r(arma::span(0, N-1));
  p = round(R*p)/R;

  // add 0 to vk
  arma::vec zero(1);
  zero.fill(0);
  vk = join_cols(zero,vk);

  // loop that select sample
  arma::vec A(N+1);
  arma::umat final(N,N);
  arma::vec tmp3(N+1);
  for(int i = 0;i < N;i++){
    tmp3.fill(cent(i));
    tmp3 = round(R*(vk-tmp3))/R;
    A = tmp3 - arma::floor(tmp3);
    final.col(i) = A(arma::span(0, N-1)) >  A(arma::span(1, N));
  }

  // remove empty last line
  final.shed_row(N-1);

  //remove duplicated value
  arma::uvec dupl = find(p < eps);
  final.shed_cols(dupl);
  final = final.t();
  arma::uvec dupl_inv = find(p >= eps);
  p = p.elem(dupl_inv);

  // be sure that pik with 1 or 0 equal 1 or 0 on each column
  arma::uvec tmp4(final.n_rows,arma::fill::ones);
  for(unsigned int k = 0;k< ones.size();k++){
    final.col(ones(k)) = tmp4;
  }
  for(unsigned int j = 0;j < zeros.size();j++){
    tmp4.fill(0.0);
    final.col(zeros(j)) = tmp4;
  }


  // return(final);
  return Rcpp::List::create(Rcpp::Named("probas")=p,
                            Rcpp::Named("samples") = final);
}')

source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/redux.R")
source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/jump.R")
source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/Orfs.R")
source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/algofastflightcube.R")
source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/fastflightcubeSPOT.R")
source("C:/Users/jauslinr/switchdrive/Sampling/Sampling/R/samplecubeSPOT.R")


set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 300
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
sampling::inclusionprobabilities(runif(N),50),
sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
dim(X)
order = 2
EPS = 1e-11

system.time(test1 <- fastflightcubeSPOT(X,pik,order = 2))
system.time(test2 <- fastflightcubeSPOT(X,pik,order = 1))
system.time(test3 <- sampling::fastflightcube(X,pik, order = 2))
system.time(test4 <- BalancedSampling::flightphase(pik,X))


system.time(s <- samplecubeSPOT(X,pik,order = 2))
system.time(s <- sampling::samplecube(X,pik,order = 2))








rm(list=ls())
user <- c('setille\\switchdrive\\__PROJETS_DE_RECHERCHE\\',
          'eustachee\\switchdrive\\',
          'jauslinr\\switchdrive\\')[3]
data      <- read.csv(file = paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\',"mitteland_square_lib.csv"), sep = ',')
df       <- data
coord    <- as.matrix(df[,1:2])
interest <- df[,3]
N        <- nrow(df)
t       <- 3
Pik     <- matrix(rep(0, t*N), ncol = t)
Pik[,1] <- rep(100/N,N)
Pik[,2] <- inclusionprobabilities(df[,4],160)
Pik[,3] <- inclusionprobabilities(df[,3],180)


S <- Orfs(Pik)





