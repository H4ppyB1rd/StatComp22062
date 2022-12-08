#include <Rcpp.h>
using namespace Rcpp;

//' @title matrix multiplication
//' @description Compute the result of matrix multiplication in cpp
//' @param matA a matrix
//' @param matB another matrix with which matA can do a matrix multiplication
//' @return The result C of C=AB
//' @examples
//' \dontrun{
//' mat1 = matrix(rnorm(15),nrow = 3,ncol = 5)
//' mat2 = matrix(rnorm(30),nrow = 5,ncol = 6)
//' c = matmultC(mat1,mat2)
//' print(c)
//' }
//' @export
//' @useDynLib StatComp22062
// [[Rcpp::export]]
NumericMatrix matmultC(NumericMatrix matA, NumericMatrix matB) {
  int m = matA.nrow(), n = matA.ncol(), p = matB.ncol();
  double stemp = 0;
  NumericMatrix matC(m, p);
  for(int i=0;i<m;++i) {
    for(int k=0;k<n;++k) {
      stemp = matA(i,k);
      for(int j=0;j<p;++j){
        matC(i,j) += stemp*matB(k,j);
      }
    }
  }
  return(matC);
}

