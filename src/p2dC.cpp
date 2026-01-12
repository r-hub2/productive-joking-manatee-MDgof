#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

//' R function order(x,y) for Rcpp
//' 
//' @param x first vector
//' @param y second vector
//' @return  a vector of integers
//' @keywords  internal
// [[Rcpp::export]]
Rcpp::IntegerVector orderC(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  int n = x.size();
  if (n != y.size()) {
    Rcpp::stop("Input vectors must have equal length!");
  }
  
  // Create a vector of indices from 0 to n-1
  std::vector<int> p(n);
  std::iota(p.begin(), p.end(), 0);
  
  // Define a custom comparator using a lambda function
  // It compares elements primarily by x, then by y for ties
  std::sort(p.begin(), p.end(), 
            [&](int i, int j){
              if (x[i] != x[j]) {
                return x[i] < x[j]; // Primary sort by x
              }
              return y[i] < y[j]; // Secondary sort by y (for ties in x)
            });
  
  // Convert the C++ index vector to an Rcpp IntegerVector (1-based indexing for R)
  Rcpp::IntegerVector res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = p[i] + 1; 
  }
  return res;
}

//' Find probabilities from cdf for discrete data
//' 
//' @param  x matrix with data
//' @param  cdf  function to find distribution function
//' @param  p (possible) arguments for cdf
//' @param  Fx (if available) already calculated values of cdf
//' @return a matrix with probabilities added
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix p2dC(Rcpp::NumericMatrix x, 
                         Rcpp::Function cdf, 
                         Rcpp::NumericVector p,
                         Rcpp::NumericVector Fx=Rcpp::NumericVector::create(-1)) {
  int k=x.nrow(),i,j;
  NumericVector x1=unique(x(_,0));
  int nx=x1.size(),ny=k/nx;
  NumericMatrix z(k, 5);
  if(Fx[0]<0) { //need to find cdf values
    Rcpp::Environment base("package:base");
    Rcpp::Function formals_r = base["formals"];
    Rcpp::List rescdf = formals_r(Rcpp::_["fun"]=cdf);
    NumericMatrix zz=z(_,Range(0, 1));
    if(rescdf.size()==1) Fx=cdf(x);
    else Fx=cdf(x, p);
  }
  IntegerVector I=orderC(x(_,0), x(_,1)); // reorder matrix
  for(i=0;i<k;++i) {
    for(j=0;j<3;++j) {
      z(i,j)=x(I[i]-1, j);
    }
    z(i, 3)=Fx[I[i]-1];
  }
  z(0, 4)=z(0, 3);
  for(i=1;i<k;++i) {
      if(i%ny==0) z(i, 4)=z(i, 3)-z(i-ny, 3);
      else {
        z(i, 4)=z(i, 3)-z(i-1, 3);
        if(i-ny>=0) z(i, 4)=z(i, 4)-z(i-ny, 3);
        if(i-ny-1>=0) z(i, 4)=z(i, 4)+z(i-ny-1, 3);
      }
  }
  return z;
}
