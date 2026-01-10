#include <Rcpp.h>
#include <cmath> // For std::sqrt, std::max
using namespace Rcpp;

//' Find test statistic for Kernel Stein Discrepancy test
//'
//' @param X data set.
//' @param scf function to find scores 
//' @param p (possible) parameters
//' @return a double (test statistic)
//' 
// [[Rcpp::export]]
double ksd(NumericMatrix X, Function scf, NumericVector p) {
  int n = X.nrow();
  int d = X.ncol();
  NumericMatrix S = scf(X, p);
  
  // Step 1: compute pairwise squared distances and (optionally) median
  double sigma2;
  std::vector<double> upper;
  upper.reserve(n*(n-1)/2);
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      double sqdist = 0.0;
      for (int k = 0; k < d; k++) {
        double diff = X(i,k) - X(j,k);
        sqdist += diff * diff;
      }
      upper.push_back(sqdist);
    }
  }
  
  std::nth_element(upper.begin(), upper.begin() + upper.size()/2, upper.end());
  double med2 = upper[upper.size()/2];
  sigma2 = 0.5 * med2;

  // Step 2: compute KSD^2 statistic
  double ksd2 = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sqdist = 0.0;
      for (int k = 0; k < d; k++) {
        double diff = X(i,k) - X(j,k);
        sqdist += diff * diff;
      }
      double kij = exp(-sqdist / (2.0 * sigma2));
      
      double term1 = 0.0, term2 = 0.0, term3 = 0.0;
      for (int k = 0; k < d; k++) {
        double diff = X(i,k) - X(j,k);
        term1 += S(i,k) * S(j,k);
        term2 += S(i,k) * (-diff) / sigma2;
        term3 += S(j,k) * ( diff) / sigma2;
      }
      double term4 = (d / sigma2 - sqdist / (sigma2 * sigma2));
      ksd2 += (term1 + term2 + term3 + term4) * kij;
    }
  }
  
  return ksd2 / ((double)n * (double)n);
}
