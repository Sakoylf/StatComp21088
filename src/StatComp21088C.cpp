#include <Rcpp.h>
using namespace Rcpp;

//' @title Rcpp function
//' @name gibbsC
//' @description an Rcpp function of gibbs sample
//' @param N sample size
//' @param a The first parameter of the distribution
//' @param b The second parameter of the distribution
//' @param n The third parameter of the distribution
//' @return generated random numbers matrix \code{mat}
//' @examples
//' \dontrun{
//' C <- gibbsC(5000,1,1,25)
//' print(C)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int a, int b, int n){
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  mat(0,0) = 0;
  mat(0,1) = 0.5;
  for(int i = 1; i < N; i++) {
    x = rbinom(1, n, mat(i-1,1))[0];
    mat(i, 0) = x;
    y = rbeta(1, x+a, n-x+b)[0];
    mat(i, 1) = y;
  }
  return(mat);
}


