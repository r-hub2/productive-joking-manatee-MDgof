#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Compute squared Euclidean distance between two vectors
double sq_dist(NumericVector x, NumericVector y) {
  int d = x.size();
  double s = 0.0;
  for(int i = 0; i < d; i++) {
    double diff = x[i] - y[i];
    s += diff * diff;
  }
  return s;
}

// Gaussian RBF kernel
double rbf_kernel(NumericVector x, NumericVector y, double sigma) {
  return exp(-sq_dist(x, y) / (2.0 * sigma * sigma));
}

// Compute median of upper triangle of pairwise squared distances
double median_sigma(NumericMatrix X) {
  int n = X.nrow();
  std::vector<double> dists;
  for(int i = 0; i < n; i++) {
    for(int j = i+1; j < n; j++) {
      dists.push_back(std::sqrt(sq_dist(X(i,_), X(j,_))));
    }
  }
  std::nth_element(dists.begin(), dists.begin() + dists.size()/2, dists.end());
  return dists[dists.size()/2];
}

// Compute MMD^2 statistic between two samples X and Y
double compute_mmd(NumericMatrix X, NumericMatrix Y, double sigma) {
  int n = X.nrow();
  int m = Y.nrow();
  
  double Kxx = 0.0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      Kxx += rbf_kernel(X(i,_), X(j,_), sigma);
  Kxx /= (n * n);
  
  double Kyy = 0.0;
  for(int i = 0; i < m; i++)
    for(int j = 0; j < m; j++)
      Kyy += rbf_kernel(Y(i,_), Y(j,_), sigma);
  Kyy /= (m * m);
  
  double Kxy = 0.0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
      Kxy += rbf_kernel(X(i,_), Y(j,_), sigma);
  Kxy /= (n * m);
  
  return Kxx + Kyy - 2.0 * Kxy;
}

//' Find test statistics for continuous data using MMD method
//' 
//' @param x A numeric matrix.
//' @param rnull routine to generate data.
//' @param p parameter for rnull
//' @return A double
// [[Rcpp::export]]
double mmd(NumericMatrix x, Function rnull, NumericVector p) {
  int n=x.nrow(), d=x.ncol();
  double sigma = median_sigma(x);
  if(sigma == 0) sigma = 1.0;
  
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  NumericMatrix Y_ref(n, d);
  if(resrnull.size()==0) Y_ref=rnull();
  else Y_ref=rnull(p);
  double stat_obs = compute_mmd(x, Y_ref, sigma);
  return stat_obs;
}
  

  
