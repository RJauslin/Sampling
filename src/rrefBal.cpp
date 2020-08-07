#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
//' @title reduced row echelon fom
//'
//'
//' @param M matrix 
//'
//' @return NULL (transform matrix)
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
void rrefBal(NumericMatrix& M){
  int lead = 0;
  int rowCount = M.nrow();
  int columnCount = M.ncol();
  double eps = 1e-11;
  int r,i,k;
  double temp;
  for(r=0; r<rowCount; r++){
    if(columnCount<=lead){return;}
    i = r;
    while( std::max(M(i,lead),-M(i,lead)) < eps ){
      M(i,lead) = 0.0;
      i = i + 1;
      if(i == rowCount){
        i = r;
        lead = lead + 1;
        if(columnCount == lead){return;}
      }
    }
    // swap rows i and r
    for(k=0;k<columnCount;k++){
      temp = M(i,k);
      M(i,k) = M(r,k);
      M(r,k) = temp;
    }
    // If M(r, lead) is not 0 divide row r by M(r, lead)
    if( M(r,lead) != 0.0 ){
      temp = M(r,lead);
      for(k=0;k<lead;k++){M(r,k) = 0.0;}
      for(k=lead;k<columnCount;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(i=0;i<rowCount;i++){
      if( i != r ){
        // Subtract M(i, lead) multiplied by row r from row i
        temp = M(i,lead);
        for(k=0;k<columnCount;k++){
          M(i,k) = M(i,k) - temp * M(r,k);
        }
      }
    }
    lead = lead + 1;
  }
  return;
}


/*** R
set.seed(1)
rm(list = ls())
N = 50
n = 30
p = 20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A <- as.matrix(X/pik)

test <- t(A[1:(p+2),])
rrefBal(test)
test

*/
