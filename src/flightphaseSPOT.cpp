#include <RcppArmadillo.h>
#include "rrefBal.h"
#include "NumericToArma.h"
#include "reduxArma.h"

using namespace Rcpp;



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
NumericVector onestepfastflightcubeSPOT(NumericVector prob, NumericMatrix Bm){
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
NumericVector flightphaseSPOT(NumericVector prob, NumericMatrix Xbal){
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
      
      arma::mat B_arma = mat_as_arma(B);
      Rcpp::List L = reduxArma(B_arma);
      NumericMatrix B_redux = L[0];
      IntegerVector ind_col = L[1];
      IntegerVector ind_row = L[2];
      
      std::cout << B_redux <<std::endl << std::endl;
      std::cout << B<<std::endl;
      
      if(howmany <= naux){
  
        
        // arma::mat B_arma = mat_as_arma(B);
        arma::mat kern = arma::null(B_arma.t());
        if(kern.empty()){
          break;
        }
      }
      
      p_small = onestepfastflightcubeSPOT(p_small,B);
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
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 50
n1 <- floor(N/3)
n2 <- floor(N/5)
n3 <- floor(N/7)
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),n1),
sampling::inclusionprobabilities(runif(N),n2),
sampling::inclusionprobabilities(runif(N),n3)),ncol = 3)
X <- PM(Pik)$PM
X <- cbind(rep(1,nrow(X)),X)
image(as(X,"sparseMatrix"))

pik <- PM(Pik)$P

dim(X)
order = 2
EPS = 1e-11

test <- flightphaseSPOT(pik,X)

test <- SamplingC::flightphase(pik,X)
test <- sampling::fastflightcube(X,pik,order = 2)
test <- BalancedSampling::flightphase(pik,X)

system.time(test <- flightphase_arma(X,pik))
system.time(test <- flightphase_arma2(X,pik))
# t(X)%*%test
A <- X/pik

t(A)%*%pik
t(A)%*%test
t(A)%*%pik # correct
t(A)%*%pikCube01



length(which(test > eps & test < (1-eps)))
length(which(test > eps & test < (1-eps)))


test <- matrix(c(1:6,rep(0,6*4)),ncol = 6,nrow = 5,byrow = T)
t <- matrix(rnorm(5*6),ncol = 6,nrow = 5)
rrefBal(t)
ukern(t)

library(MASS)
Null(t(test))

rm(list = ls())
N = 50
n = 30
p = 2
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
flightphaseSPOT(pik,X)
*/