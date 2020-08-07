#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Arma matrix to NumericMatrix
//'
//' @param x a Matrix
//'
//' @return same matrix but with NumericMatrix type.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix mat_as_Numeric (arma::mat x) {
  Rcpp::NumericMatrix y = Rcpp::wrap(x) ;
  return(y) ;
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title NumericMatrix matrix to Arma
//'
//' @param x a Matrix 
//'
//' @return same matrix but with arma type.
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::mat mat_as_arma(Rcpp::NumericMatrix x) {
  arma::mat y = Rcpp::as<arma::mat>(x);
  return(y) ;
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title svd in arma
//'
//' @param x Matrix X
//'
//' @return list of element of svd
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List svdArma(arma::mat x) {
  arma::mat U ;
  arma::vec s ;
  arma::mat V ;
  svd(U, s, V, x,"dc") ;
  Rcpp::List out ;
  out["U"] = U ;
  out["s"] = s ;
  out["V"] = V ;
  return(out) ;
}

/*** R

X <- matrix(rnorm(10000),ncol = 100,nrow = 100)
system.time(mat_as_arma(X))
system.time(mat_as_Numeric(X))



X <- matrix(rnorm(10000),ncol = 100,nrow = 100)
system.time(test <- svdArma(X))
system.time(test2 <- svd(X))

library(microbenchmark)
set.seed(42)
X <- matrix(rnorm(25e4), 5e2, 5e2)
microbenchmark(svd(X), svdArma(X))


*/
