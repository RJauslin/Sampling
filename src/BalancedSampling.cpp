#include <RcppArmadillo.h>
#include "NumericToArma.h"
#include "rrefBal.h"
using namespace Rcpp;

//**********************************************
// Author: Anton Grafström
// Last edit: 2014-05-05
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
bool all0(NumericVector x) {
  double eps = 1e-12;
  return is_true(all((x < eps) & (x > -eps)));
}

/*** R
N <- 10
all0(c(rep(0,N),1))
*/


// [[Rcpp::depends(RcppArmadillo)]]
//' @title onestepfastflightcube
//'
//' @description
//' 
//' one step of the fast flight cube. Direct implementation of the package BalancedSampling from Anton Grafström
//'
//' @param prob vector of inclusion probabilities
//' @param Bm matrix of size (p x p+1) (transpose from the original)
//'
//' @return updated inclusion probabilities
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
  
  // u = ukern(Bm);
  
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
//' @title Flight phase of the cube method
//'
//' @description
//' 
//' Direct implementation of the function \code{\link[BalancedSampling:flightphase]{flightphase}}. 
//' 
//' In this function there exist a break condition when the number of not updated
//' inclusion probabilities reach p (number of auxiliary variables).
//' 
//' It is possible that in fact the kernel of the matrix is not empty and then the algorithm could
//' continue without loosing any information.
//' 
//' Look at \code{\link{ffphase}} for a modified version of this function that continue until we have an empty kernel.
//' (and still have a good time-consuming)
//'
//' @param prob vector of length N with inclusion probabilities
//' @param Xbal matrix of balancing auxiliary variables of N rows and q columns
//'
//' @return Returns a vector of length N with new probabilities, where at most q are non-integer.
//'
//' @export
// [[Rcpp::export]]
NumericVector flightphase(NumericVector prob, NumericMatrix Xbal){
  int N = prob.size();
  int naux = Xbal.ncol();
  
  IntegerVector index(N);
  NumericVector p(N);
  int i,j,howmany;
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
    if(howmany <= naux){done=N; break;} //WHY ?!?!?!
    
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
        
        // check singular value --> if all are greater than 0 (eps)
        // then the matrix has a empty kernel -> IN FACT LOWER (certainly due to the all function)
        // arma::vec s = arma::svd(B_arma);
        // std::cout << s << std::endl;
        // if(arma::all(s > eps)){
        // break;
        // }
        
        
        // NumericVector u = ukern(B);
        // if(sum(u) < eps){
        //   break;
        // }
        
        // arma::mat B_arma = mat_as_arma(B);
        // arma::mat kern = arma::null(B_arma.t());
        // std::cout << B << std::endl;
        // if(kern.empty()){
        //   emp = true;
        //   // break;
        // }
      }
      // std::cout << B << std::endl;
      
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
  
  
  // arma::mat B_arma = mat_as_arma(B);
  // arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
  // 
  
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