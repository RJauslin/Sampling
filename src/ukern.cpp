#include <RcppArmadillo.h>
#include "rrefBal.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title kernel vector
//' @param Bm matrix
//' @return a null vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @export
// [[Rcpp::export]]
NumericVector ukern(NumericMatrix Bm){
  int ncol = Bm.ncol();
  int nrow = Bm.nrow();
  
  NumericVector u(ncol,0.0);
  IntegerVector uset(ncol,0);
  double eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  rrefBal(Bm);
  
  // std::cout << Bm << std::endl;
  for(int i = (nrow-1); i >=0; i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(int j=0;j < ncol; j++){
      if(Bm(i,j)==0.0){
        lead++;
      }else{
        break;
      }
    }
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(int j = lead + 1;j < ncol; j++){
        if( uset[j] == 0 ){
          uset[j] = 1;
          free *= -1.0;
          u[j] = free;
        }
        v -= u[j]*Bm(i,j);
      }
      u[lead] = v/Bm(i,lead);
      uset[lead] = 1;
    }
  }
  // unset u[i] are free and are set to 1 or -1, can only exist at beginning
  for(int i=0;i<ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{break;}
  }
  
  return u;
}
/*** R
# set.seed(1)
rm(list = ls())
N = 50
n = 30
p = 20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A <- as.matrix(X/pik)

test <- t(A[1:(p+2),])
rankMatrix(test)
dim(test)
u <- ukern(test)

test
u
*/


