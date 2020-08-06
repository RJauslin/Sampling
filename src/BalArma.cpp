#include <RcppArmadillo.h>
#include "rrefArma.h"

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
arma::vec onestepfastflightcubeArma(arma::vec prob, arma::mat Bm){
  int ncol = Bm.n_cols;
  int nrow = Bm.n_rows;
  
  arma::vec u(ncol,arma::fill::zeros);
  arma::uvec uset(ncol,arma::fill::zeros);

  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  rrefArma(Bm);
  
  // std::cout << Bm << std::endl;
  for(int i = (nrow-1);i >= 0; i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(int j = 0; j < ncol; j++){
      if(Bm(i,j)==0.0){
        lead++;
      }else{
        break;
      }
    }
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(int j = lead+1;j < ncol;j++){
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
  for(int i = 0;i < ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{
      break;
    }
  }
  
  // u = ukern(Bm);
  
  // find lambda1 and lambda2
  for(int i = 0;i < ncol;i++){
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
  for(int i = 0;i < ncol;i++){
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
arma::vec flightphaseArma(arma::vec prob, arma::mat Xbal){
  int N = prob.size();
  int naux = Xbal.n_cols;
  
  arma::uvec index(N);
  arma::vec p(N);
  
  int i,j,k,howmany;
  for(int i = 0;i < N;i++){
    index[i]=i;
    p[i]=prob[i];
  }
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
  for(int i = done;i < N;i++){
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
    howmany = std::min(naux + 1,N-done);
    
    // stop if there are less than naux units left
    // if(howmany <= naux){
    //   done=N;
    //   break;
    // } //WHY ?!?!?!
    
    if( howmany > 1 ){
      arma::vec p_small(howmany);
      arma::vec dists(howmany); dists = 1e+20;
      arma::vec index_small(howmany);
      // arma::mat B(howmany-1,howmany);
      arma::mat B(naux,howmany);
      
      
      for(int i = 0;i < howmany; i++){
        index_small[i] = index[done+i];
        for(int j = 0;j < naux;j++){
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }
      
      if(howmany < naux + 1){
        arma::mat kern = arma::null(B);
        if(kern.empty()){
          break;  
        }
      }
      
        //   emp = true;
        //   // break;
        // }
      // std::cout << index << std::endl;
      // std::cout << B << std::endl;
      // if(howmany <= naux){
        
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
      // }
      // std::cout << B << std::endl;
      
      p_small = onestepfastflightcubeArma(p_small,B);
      
      
      // update prob
      for(int i = 0;i < howmany;i++){
        p[index_small[i]] = p_small[i];
      }
      // update done and index
      howlong = done + howmany;
      for(int i = done;i < howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
      if(runif(1)[0] < p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }
  }
  
  
  // round
  for(int i = 0;i < N;i++){
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
flightphaseArma(pik,X)





rm(list = ls())
N = 5000
n = 800
p = 40
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
system.time(test1 <- BalancedSampling::flightphase(pik,X))
system.time(test2 <- flightphaseArma(pik,X))


A <- X/pik

t(A)%*%pik
t(A)%*%test1
t(A)%*%test2
t(A)%*%pik # correc









  
*/