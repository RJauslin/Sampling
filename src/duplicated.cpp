#include <RcppArmadillo.h>

using namespace Rcpp;


//' @title all test condition in Rcpp
//'
//' @param x vector of boolean
//'
//' @return true or false
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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


//' @title colSums in Rcpp
//'
//' @param x Matrix X
//'
//' @return colSums
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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


//' @title rowSums in Rcpp
//'
//' @param x Matrix X
//'
//' @return rowSums
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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



//' @title Check equal
//'
//' @param x vector 
//' @param y vector
//'
//' @return return true if vector equal to each other
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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



//' @title duplicated columns
//'
//' @param B Matrix B
//'
//' @return return vector of boolean
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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


//' @title duplicated rows
//'
//' @param B Matrix B
//'
//' @return vector of boolean
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
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


