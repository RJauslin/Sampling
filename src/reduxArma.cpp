#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Reduction of the matrix
//'
//' @description 
//' 
//' This function reduces the size of the matrix by removing alternatively columns and rows that have sum equal to 0.
//' 
//' In case where the number of auxiliary varibale is great (p very large), even if we use the fast implementation proposed by
//' (Chauvet and Tillé 2005) the problem is time consuming. If we have the chance that the matrix is strongly sparse,
//' we can then use the function to reduce the size of the matrix B by using this method. 
//' 
//' If the matrix is dense or every column have sum greater than 0, then nothing is changed.
//'
//' @param B a matrix of size (p+1 x p) sub-matrix of auxiliary matrix.
//'
//' @return a list
//' 
//' @references 
//' Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53–62. 
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List reduxArma(arma::mat B) {
  
  double eps = 1e-12;
  
  arma::mat B_out = B;
  arma::rowvec sums = sum(B,0);
  arma::colvec sums_row = sum(B,1);
  
  arma::uvec ind_col = arma::regspace<arma::uvec>(0,1,B_out.n_cols-1);
  arma::uvec ind_row = arma::regspace<arma::uvec>(0,1,B_out.n_rows-1);
  
  int step = 1;
  while(any(sums < eps && sums > -eps)){ // loop while we find some colSums equal to 0
    
    // std::cout << step << std::endl << std::endl;

    
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
    
   
    
    
    // std::cout << B_out.n_cols << std::endl;
    // std::cout << B_out.n_rows << std::endl;
    
    
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
                              Rcpp::Named("ind_col") = ind_col ,
                              Rcpp::Named("ind_row") = ind_row );
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