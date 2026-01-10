#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistic of Fasano–Franceschini test
//' @param dta data matrix
//' @return  a test statistic
//' @keywords internal
// [[Rcpp::export]]
double FF(NumericMatrix dta) {
  NumericVector x=dta(_,0), y=dta(_,1);
  int n = dta.nrow();
  double Dmax = 0.0;
  
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    double yi = y[i];
    
    int c_ll = 0, c_ul = 0, c_lr = 0, c_ur = 0;
    
    for (int j = 0; j < n; ++j) {
      if (x[j] <= xi && y[j] <= yi) c_ll++;
      else if (x[j] <= xi && y[j] > yi) c_ul++;
      else if (x[j] > xi && y[j] <= yi) c_lr++;
      else c_ur++;
    }
    
    double p_ll = (double)c_ll / n;
    double p_ul = (double)c_ul / n;
    double p_lr = (double)c_lr / n;
    double p_ur = (double)c_ur / n;
    
    double e_ll = xi * yi;
    double e_ul = xi * (1.0 - yi);
    double e_lr = (1.0 - xi) * yi;
    double e_ur = (1.0 - xi) * (1.0 - yi);
    
    Dmax = std::max(Dmax, std::abs(p_ll - e_ll));
    Dmax = std::max(Dmax, std::abs(p_ul - e_ul));
    Dmax = std::max(Dmax, std::abs(p_lr - e_lr));
    Dmax = std::max(Dmax, std::abs(p_ur - e_ur));
  }
  
  return Dmax;
}

  
