#include <RcppArmadillo.h>
#include "NumericToArma.h"
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-05-05
// Licence: GPL (>=2)
//**********************************************

// "import" print for error checking
Function print("print");


// [[Rcpp::export]]
bool all0(NumericVector x) {
  double eps = 1e-12;
  return is_true(all(x < eps & x > -eps));
}

/*** R
N <- 10
all0(c(rep(0,N),1))
*/


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
NumericVector onestepfastflightcube(NumericVector prob, NumericMatrix Bm){
  int ncol = Bm.ncol();
  int nrow = Bm.nrow();
  int i, j;
  NumericVector u(ncol,0.0);
  IntegerVector uset(ncol,0);
  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  rrefBal(Bm);
  
  // std::cout << Bm << std::endl;
  for(i=(nrow-1);i>=0;i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(j=0;j<ncol;j++){if(Bm(i,j)==0.0){lead++;}else{break;}}
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(j=lead+1;j<ncol;j++){
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
  for(i=0;i<ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{break;}
  }
  
  // find lambda1 and lambda2
  for(i=0;i<ncol;i++){
    if(u[i]>0){
      la1 = std::min(la1,(1-prob[i])/u[i]);
      la2 = std::min(la2,prob[i]/u[i]);
    }
    if(u[i]<0){
      la1 = std::min(la1,-prob[i]/u[i]);
      la2 = std::min(la2,(prob[i]-1)/u[i]);
    }
  }
  // random choice of p+lambda1*u and p-lambda2*u
  if(runif(1)[0]<la2/(la1+la2)){
    la = la1;
  }else{
    la = -la2;
  }
  // update prob
  for(i=0;i<ncol;i++){
    prob[i] = prob[i] + la * u[i];
    if(prob[i] < eps){ prob[i] = 0; }
    if(prob[i] > 1-eps){ prob[i] = 1; }
  }
  return prob;
}


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
NumericVector flightphase(NumericVector prob, NumericMatrix Xbal){
  int N = prob.size();
  int naux = Xbal.ncol();
  
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}
  double eps = 1e-12;
  int done = 0, tempInt, howlong;
  // randomize order of index list
  // NumericVector rnd = runif(N);
  // for(i=0;i<N;i++){
  //   k = i + floor(rnd[i] * (N-i));
  //   tempInt = index[i];
  //   index[i] = index[k];
  //   index[k] = tempInt;
  // }
  
  // put finished units at beginning of list
  for(i=done;i<N;i++){
    if( p[index[i]]<eps || p[index[i]]>1-eps ){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  }
  
  // remaining are index from done to N-1
  while( done < N ){
  // while( B.ncol() > 1){
    
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    
    // stop if there are less than naux units left
    // if(howmany <= naux){done=N; break;} //WHY ?!?!?!
    
    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      NumericMatrix B(howmany-1,howmany);
      for(i=0;i<howmany;i++){
        index_small[i] = index[done+i];
        for(j=0;j<howmany-1;j++){
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }
      // std::cout << index << std::endl;
      // std::cout << B << std::endl;
      if(howmany <= naux){
        
        // NumericVector u = ukern(B);
        // std::cout << u << std::endl;
        // if(all0(u) == true){
        //   break;
        // }
        
        
        // std::cout << "B.ncol smaller" << std::endl;
        arma::mat B_arma = mat_as_arma(B);
        arma::mat kern = arma::null(B_arma.t());
        if(kern.empty()){
          break;
        }
      }
      
      p_small = onestepfastflightcube(p_small,B);
      // update prob
      for(i=0;i<howmany;i++){
        p[index_small[i]] = p_small[i];
      }
      // update done and index
      howlong = done + howmany;
      for(i=done;i<howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
      if(runif(1)[0]<p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }
  }
  // round
  for(i=0;i<N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  
  return p;
}

/*** R
rm(list = ls())
N = 50
n = 30
p = 2
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
flightphase(pik,X)
*/