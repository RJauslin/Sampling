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
  while(any(sums > eps || sums < -eps)){
    std::cout << step << std::endl;
    // intialize column an row index
    ind_col = ind_col.elem(arma::find(sums > eps || sums < -eps));
    
    
    // extract B_out
    B_out = B_out.cols(ind_col);
    
    
    // calculate unique rows_sums and unique rows
    sums_row = sum(B_out,1);
    ind_row = ind_row.elem(arma::find(sums_row > eps || sums_row < -eps)); // find > 0 or < 0
    B_out = B_out.rows(ind_row); // update B_out
    sums_row = sum(B_out,1); // recompute sums_row
    ind_row = ind_row.elem(arma::find(unique(sums_row)));// find unique value
    
    
    // std::cout << ind_row << std::endl;
    
    if(ind_row.size() > (B_out.n_cols + 1)){
      arma::uvec f = arma::regspace<arma::uvec>(0,1,B_out.n_cols);
      B_out = B_out.rows(ind_row.elem(f));
    }else{
      if(ind_row.size() == B_out.n_cols){
        arma::uvec c = arma::regspace<arma::uvec>(0,1,B_out.n_cols-2);
        arma::uvec r = arma::regspace<arma::uvec>(0,1,B_out.n_cols-1);
        
        ind_row = ind_row.elem(r);
        ind_col = ind_col.elem(c);
        
        B_out = B_out(r,c);
      }
      break;
    }
    
    sums = sum(B_out,0);
    step = step + 1;
  }

  
  return Rcpp::List::create(Rcpp::Named("B") = B_out,
                            Rcpp::Named("ind_col") = ind_col +1,
                            Rcpp::Named("ind_row") = ind_row +1);
  // return(B_out);
}


/*** R
rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(100,100,density = 0.001))
image(as(B,"sparseMatrix"))
reduxArma(B)
reduxB(B)




rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(200,200,density = 0.001))
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
system.time(test2 <- reduxB(B))

*/


//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
bool all_sug(LogicalVector x) {
  return is_true(all(x == TRUE));
}



/*** R
x <- c(3, 9, 0, 2, 7, 5, 6)
y <- c(0, 0, 0, 0, 0, 0, 0)
all_sug(y == 0)
*/


//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
NumericVector colSumsRcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += x(i, j);
    }
    ans[j] = sum;
  }
  return ans;
}
/*** R
set.seed(1)
X <- as.matrix(rsparsematrix(100,100,density = 0.01))
system.time(test <- colSums(X))
system.time(test1 <- colSumsRcpp(X))

*/


//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
NumericVector rowSumsRcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}
/*** R
set.seed(1)
X <- as.matrix(rsparsematrix(4000,4000,density = 0.01))
system.time(colSums(X))
system.time(colSumsRcpp(X))

*/

//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
NumericVector reduxRcpp(NumericMatrix B) {
  double eps = 1e-12;
  int N = B.nrow();
  int p = B.ncol();
  NumericVector sums = colSums(B);
  NumericVector sums_row = rowSums(B);
  bool check = false; 
  
  
  IntegerVector ind_col(p);
  IntegerVector ind_row(N);
  for(int i=0; i < p; i++){
    ind_col[i]=i;
  }
  for(int i=0; i < N; i++){
    ind_row[i]=i;
  }
  
  
  while(check == false){
    
    check = is_true(any(sums < eps | sums > -eps));
  }
  
  return(sums_row);
}


/*** R
set.seed(1)
X <- as.matrix(rsparsematrix(100,100,density = 0.01))
reduxRcpp(X)


*/





//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
NumericVector duplicatedMatrixRcpp(NumericVector x) {
  return Rcpp::unique(x);
}



/*** R
x <- c(3, 0, 0, 2, 7, 5, 6)
unique(x)
*/


//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
bool check_equal(NumericVector x, NumericVector y) {
  bool out = false;
  if (all_sug(x == y)) {
    out = true;
  } 
  return(out);
}



/*** R
x <- c(3, 9, 0, 2, 7, 5, 6)
y <- c(0, 0, 0, 0, 0, 0, 0)
# y <- c(3, 9, 0, 2, 7, 5, 6)
check_equal(x,y)
*/



//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
LogicalVector duplicatedCols(NumericMatrix B) {
  int N = B.nrow();
  int p = B.ncol();
  bool check_tmp = false;
  LogicalVector check(p,false);
  for(int i = 0; i< p; i ++){
    for(int j = i+1; j < p ;j++ ){
      check_tmp = check_equal(B.column(i),B.column(j));
      if(check_tmp == true){
        check[i] = true;
        if(j == p-1){
          check[j]= true;
        }
        check_tmp = false;
      }
    }
  }
  return(check);
}
/*** R
B <- matrix(c(1:6,1:6,1:6),ncol = 3)
duplicatedCols(B)
B <- matrix(c(1:6,1:6,1:6,7:12,1:6),ncol = 5)
duplicatedCols(B)




set.seed(1)
X <- as.matrix(rsparsematrix(1000,1000,density = 0.01))
system.time(duplicatedCols(X))
system.time(duplicated(X,MARGIN = 2))

*/


//' @title Arma matrix to NumericMatrix
//'
//' @param x Matrix X
//'
//' @return same matrix
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//
//'
//' @export
// [[Rcpp::export]]
LogicalVector duplicatedRows(NumericMatrix B) {
  int N = B.nrow();
  bool check_tmp = false;
  LogicalVector check(N,false);
  for(int i = 0; i< N; i ++){
    for(int j = i+1; j < N ;j++ ){
      check_tmp = check_equal(B.row(i),B.row(j));
      if(check_tmp == true){
        check[i] = true;
        if(j == N-1){
          check[j]= true;
        }
        check_tmp = false;
      }
    }
  }
  return(check);
}
/*** R
B <- matrix(c(1:6,1:6,1:6),nrow = 3,byrow = TRUE)
duplicatedRows(B)
B <- matrix(c(1:6,1:6,1:6,7:12,1:6),nrow = 5,byrow = TRUE)
duplicatedRows(B)
duplicated(B,MARGIN = 1)







*/






// [[Rcpp::export]]
bool is_duplicate_row(R_xlen_t r, LogicalMatrix x) {
  R_xlen_t i = 0, nr = x.nrow();
  const LogicalMatrix::Row& y = x.row(r);
  
  for (; i < r; i++) {
    if (is_true(all(y == x.row(i)))) {
      return true;
    }
  }
  for (i = r + 1; i < nr; i++) {
    if (is_true(all(y == x.row(i)))) {
      return true;
    }
  }
  
  return false;
}





// Rcpp::NumericMatrix removeDuplicatedRows(NumericMatrix B) {
//   //initialization
//   return
//   
// }
// 

// Rcpp::List reduxRcpp(NumericMatrix B) {
//   //initialization
//   double eps = 1e-12;
//   sums <- colSums(B)
//   sums_row <- rowSums(B)
//   B_out <- B
//   ind_col <- seq(1,ncol(B_out),1)
//   ind_row <- seq(1,nrow(B_out),1)
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//