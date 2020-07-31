#ifndef NumericToArma_H
#define NumericToArma_H

#include <RcppArmadillo.h>
Rcpp::NumericMatrix mat_as_Numeric (arma::mat x);
arma::mat mat_as_arma(Rcpp::NumericMatrix x);
Rcpp::List svdArma(arma::mat x);

#endif
