#include <RcppArmadillo.h>
  
  // [[Rcpp::depends(RcppArmadillo)]]
  //' @title title
  //'
  //' @description
  //' description
  //'
  //'
  //' @param x x
  //'
  //' @details
  //'
  //' details
  //'
  //' @return a vector
  //'
  //'
  //' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
  //'
  //' @seealso
  //' func
  //'
  //' @examples
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