#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix gibbsC(double mu_x,double mu_y, double sig_x,double sig_y, double rho, double init_x, double init_y,int N) {
  NumericMatrix mat(N, 2);
  double s_x = sqrt(1-rho*rho)*sig_x;
  double s_y = sqrt(1-rho*rho)*sig_y;
  double x = init_x, y = init_y;
  double m1,m2;
  for(int i = 0; i < N; i++) {
    m1 = mu_x+rho*(y-mu_y)*sig_x/sig_y;
    x = rnorm(1,m1,s_x)[0];
    m2 = mu_y+rho*(x-mu_x)*sig_y/sig_x;
    y = rnorm(1,m2,s_y)[0];
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
