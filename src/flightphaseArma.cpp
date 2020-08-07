#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param B matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return updated pik
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec onestep(arma::mat B,arma::vec pik,double EPS=0.0000001){
  
  arma::mat kern = arma::null(B);
  arma::uword N = pik.size();
  arma::vec u(N);
  u = kern.col(0);
  // int ncol = kern.n_cols;
  double l1 = 1e+200;
  double l2 = 1e+200;
  double l = 1e-9;
  
  for(arma::uword k = 0; k < N; k++){
    if(u[k]> 0){
      l1 = std::min(l1,(1.0 - pik[k])/u[k]);
      l2 = std::min(l2,pik[k]/u[k]);
    }
    if(u[k]< 0){
      l1 = std::min(l1,-pik[k]/u[k]);
      l2 = std::min(l2,(pik[k]-1.0)/u[k]);
    }
  }
  if(Rcpp::runif(1)[0]<l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  for(arma::uword k = 0; k < N; k++){
    pik[k] = pik[k] + l*u[k];
    if(pik[k] < EPS){
      pik[k] = 0;
    }
    if(pik[k] > (1-EPS)){
      pik[k] = 1;
    }
  }
  return(pik);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//'
//'
//' @param X matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return a sample
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec flightphase_arma2(arma::mat X,arma::vec pik,double EPS=0.0000001){
  
  /*
   * Initialization of constant of the function
   */
  unsigned int J = X.n_cols;
  unsigned int N = pik.size();
  arma::uvec index(N,arma::fill::zeros);
  int work = 0;
  arma::uword changed = 0;
  arma::vec p(N,arma::fill::zeros);
  for(arma::uword k=0;k < N;k++){
    index[k]=k;
    p[k]=pik[k];
  }
  
  /*
   * Put done unit at the begining in index -> same as Grafström
   */
  int tmp = 0;
  for(arma::uword l = 0; l < N; l++){
    if( pik[index[l]]<EPS || pik[index[l]]>1-EPS ){
      tmp = index[work];
      index[work] = index[l];
      index[l] = tmp;
      work = work + 1;
    }
  }
  
  /*
   * Initialisation of Bsize
   */
  arma::uword Bsize = std::min(J+1,N-work);
  
  /*
   * Main LOOP done until break by null space empty
   */
  while(Bsize > 1){
    
    /*
     * sub-pik and sub-matrix B initialization
     */
    arma::vec index_work(Bsize,arma::fill::zeros);
    arma::vec pik_work(Bsize,arma::fill::zeros);
    // arma::mat B(Bsize-1,Bsize,arma::fill::zeros);
    arma::mat B(J,Bsize,arma::fill::zeros);
    
    
    /*
     * fill B with X/pik and get the right index_work
     */
    for(arma::uword i = 0;i < Bsize; i++){
      index_work[i] = index[work+i];
      // for(arma::uword j = 0; j < (Bsize-1); j++){
      for(arma::uword j = 0; j < J; j++){
        B(j,i) = X(index_work[i],j)/p[index_work[i]];
      }
      pik_work[i] = pik[index_work[i]];
    }
    
    /*
     * check if the Bsize is small that the cluster size and
     * stop if null space empty.
     */
    if(Bsize < (J+1)){
      arma::mat kern = arma::null(B);
      if(kern.empty()){
        break;
      }
    }
    // std::cout << B.n_cols << std::endl;
    // std::cout << B.n_rows << std::endl << std::endl;
    
    /*
     * onestep flightphase with the sub-pik and sub-matrix B
     */
    pik_work =  onestep(B,pik_work);
    
    /*
     * update the pik at the right position -> index_work
     */
    for(arma::uword i = 0;i < Bsize;i++){
      pik[index_work[i]] = pik_work[i];
    }
    
    // update done and index
    changed = work + Bsize;
    for(arma::uword l = work; l < changed; l++){
      if( pik[index[l]]<EPS || pik[index[l]]>1-EPS ){
        tmp = index[work];
        index[work] = index[l];
        index[l] = tmp;
        work = work + 1;
      }
    }
    
    /*
     * redefine the Bsize
     */
    Bsize = std::min(J+1,N-work);
  }
  
  /*
   * round value that are 1e-3000 or 0.9999999 by a "real" 0 or 1.
   */
  for(arma::uword i = 0;i < N;i++){
    if( pik[index[i]] > 1-EPS  ){
      pik[index[i]] = 1;
    }
    if( pik[index[i]] < EPS  ){
      pik[index[i]] = 0;
    }
  }
  return(pik);
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//'
//'
//' @param X matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return a sample
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec flightphase_arma(arma::mat X,arma::vec pik,double EPS=0.0000001){
  
  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B
  
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat B = (A.rows(i)).t(); // extract B of A
  
  while(i.size() > 0){
    // std::cout << i.size() << std::endl;
    pik.elem(i) =  onestep(B,pik.elem(i));
    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    B = (A.rows(i)).t();
    // std::cout << B.n_cols << std::endl;
    // std::cout << B.n_rows << std::endl << std::endl;
    if(i.size() < (J+1)){
      arma::mat kern = arma::null(B);
      if(kern.empty()){
        break;
      }
    }
  }
  return(pik);
}



/*** R
library(sampling)
rm(list = ls())
N = 5000
n = 200
p = 100
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
system.time(test <- flightphase_arma(X,pik))
t(A)%*%test
t(A)%*%pik
system.time(pikstar <- fastflightcubeSPOT(X, pik))
t(A)%*%pikstar
t(A)%*%pik

system.time(test <- IneqCube::flightphase(pik,X))
system.time(test <- BalancedSampling::flightphase(pik,X))
system.time(test <- fast.flight.cube(X,pik))







rm(list = ls())
 set.seed(1)
 eps <- 1e-13
 library(Matrix)
 N <- 250
 Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
 sampling::inclusionprobabilities(runif(N),50),
 sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
 X <- PM(Pik)$PM
 pik <- PM(Pik)$P
 dim(X)
 order = 2
 method = 1
 comment = TRUE
 EPS = 1e-11

 system.time(test <- flightphase_arma(X,pik))
 # t(X)%*%test
 # t(X)%*%pik
 system.time(pikstar <- fastflightcubeSPOT(X, pik))
 # t(X)%*%pikstar
 # t(X)%*%pik
 system.time(pikstar2 <- sampling::fastflightcube(X, pik))
 # t(X)%*%pikstar2
 # t(X)%*%pik
 system.time(test <- BalancedSampling::flightphase(pik,X))
 system.time(test <- fast.flight.cube(X,pik))
 t(X)%*%test
 t(X)%*%pik
 
*/