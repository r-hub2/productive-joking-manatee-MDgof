#include <Rcpp.h>
using namespace Rcpp;

//' Find gradient of log(f) 
//' @param x point of evaluation
//' @param f function 
//' @return  a gradient vector
//' @keywords internal
// [[Rcpp::export]]
NumericVector grad_vec(NumericVector x, Function f) {
  double h=0.00001;
  int i, d=x.size();
  NumericVector grad(d);
  NumericVector tmp0, tmp1, tmp2;
  tmp0=f(x);
  for(i=0;i<d;++i) {
    if(tmp0[0]<h) {
       grad[i]=0.0;
    }
    else {
      NumericVector z1 = clone(x), z2 = clone(x);
      z1[i]=z1[i]+h;
      z2[i]=z2[i]-h;
      tmp1=f(z1);
      tmp2=f(z2);
      grad[i]=(tmp1[0]-tmp2[0])/2.0/h/tmp0[0];
    }  
  }
  return grad;
}

//' Find gradient of log(f) for a matrix of points
//' @param x point of evaluation
//' @param f function 
//' @return  a matrix of gradient vectors
// [[Rcpp::export]]
NumericMatrix grad_mat(NumericMatrix x, Function f) {
   int i, n=x.nrow(), d=x.ncol();
   NumericMatrix grad(n, d);
   for(i=0;i<n;++i) {
      grad(i,_)=grad_vec(x(i,_), f);
   }
   return grad;
 }


