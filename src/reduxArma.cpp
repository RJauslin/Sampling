#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Column sums for sparseMatrix
//'
//' @description
//' Form column sums for sparseMatrix.
//'
//' @param x A sparse matrix, i.e., inheriting from \code{\link[Matrix]{sparseMatrix}}.
//'
//' @details
//' This function is designed to be used for internal \code{RcppArmadillo} functions. Nevertheless it could be applied in R.
//' It loops on the non-zero entries of the \code{\link[Matrix]{sparseMatrix}}. For general uses, the function
//' \code{\link[Matrix]{colSums}} should be prefered.
//'
//' @return column sums of x.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @seealso
//' \code{\link[Matrix]{colSums}}, \code{\link[Matrix]{rowSums}}.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List reduxArma(arma::mat B) {
  
  double eps = 1e-12;
  int N = B.n_rows;
  int p = B.n_cols;
  
  arma::mat B_out = B;
  arma::rowvec sums = sum(B,0);
  arma::colvec sums_row = sum(B,1);
  
  arma::uvec ind_col = arma::regspace<arma::uvec>(0,1,B_out.n_cols-1);
  arma::uvec ind_row = arma::regspace<arma::uvec>(0,1,B_out.n_rows-1);
  
  int step = 1;
  while(any(sums < eps && sums > -eps)){ // loop while we find some colSums equal to 0
  // while(step < 4){ // loop while we find some colSums equal to 0
    
    std::cout << step << std::endl << std::endl;
    

    
    // extract the column that have sums greater than 0
    arma::uvec coltmp = arma::find(sums > eps || sums < -eps);
    if(coltmp.size() <= 1){
      break;
    }
    B_out = B_out.cols(coltmp);
    ind_col = ind_col.elem(coltmp); // keep right index of B
    
    // calculate rowSums that have sums greater than 0
    sums_row = sum(B_out,1);
    arma::uvec rowtmp = arma::find(sums_row > eps || sums_row < -eps);
    B_out = B_out.rows(rowtmp);
    ind_row = ind_row.elem(rowtmp); // keep rignt index of B
    
    // recompute rowSums
    sums_row = arma::sum(B_out,1); 
    
   
    
    
    std::cout << B_out.n_cols << std::endl;
    std::cout << B_out.n_rows << std::endl;
    
    
    if(B_out.n_rows >= (B_out.n_cols + 1)){
      arma::uvec f = arma::regspace<arma::uvec>(0,1,B_out.n_cols); // from 0 to B_out.n_cols so B_out.n_cols + 1 element
      ind_row = ind_row.elem(f);
      B_out = B_out.rows(f);
    }else{
      // std::cout << ind_row.size() << std::endl;
      // std::cout << B_out.n_cols << std::endl;
      if(B_out.n_rows == B_out.n_cols){
        arma::uvec c = arma::regspace<arma::uvec>(0,1,B_out.n_cols-2);
        arma::uvec r = arma::regspace<arma::uvec>(0,1,B_out.n_cols-1);
        
        ind_row = ind_row.elem(r);
        ind_col = ind_col.elem(c);
        
        B_out = B_out(r,c);
      }else{
        arma::uvec r = arma::regspace<arma::uvec>(0,1,B_out.n_rows - 1);
        arma::uvec c = arma::regspace<arma::uvec>(0,1,B_out.n_rows - 2);
        
        ind_row = ind_row.elem(r);
        ind_col = ind_col.elem(c);
        
        B_out = B_out(r,c);
      }
    }
    
    sums = arma::sum(B_out,0); // update colSums
    step = step + 1;
    if(step > 100){
      break;
    }
  }
    
    return Rcpp::List::create(Rcpp::Named("B") = B_out,
                              Rcpp::Named("ind_col") = ind_col +1,
                              Rcpp::Named("ind_row") = ind_row +1);
  }
  
  
/***R
  
rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 50
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),5),
sampling::inclusionprobabilities(runif(N),10),
sampling::inclusionprobabilities(runif(N),15)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
dim(X)
order = 2
EPS = 1e-11


p <- ncol(X)
A <-  X/pik
B <- A[1:(p + 1), ]

tmp <- reduxArma(B)
tmp2 <- reduxB(B)
B_redux <- tmp$B
B_redux[1:10,1:20]
B[tmp$ind_row[1:10],tmp$ind_col[1:20]]




rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(200,200,density = 0.01))
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
dim(test1$B)
system.time(test2 <- reduxB(B))
dim(test2$B)



rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(200,200,density = 0.001))
B <- cbind(rep(1,200),B)
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
test1
dim(test1$B)
system.time(test2 <- reduxB(B))
test2
dim(test2$B)







rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 50
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),5),
                sampling::inclusionprobabilities(runif(N),10),
                sampling::inclusionprobabilities(runif(N),15)),ncol = 3)
X <- PM(Pik)$PM
image(as(X,"sparseMatrix"))
pik <- PM(Pik)$P
dim(X)
order = 2
EPS = 1e-11

system.time(test1 <- reduxArma(X))
system.time(test2 <- reduxB(X))




rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(4000,3000,density = 0.0001))
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
dim(test1$B)
system.time(test2 <- reduxB(B))
dim(test2$B)
all(as.vector(test1$ind_row) == as.vector(test2$ind_row))
*/