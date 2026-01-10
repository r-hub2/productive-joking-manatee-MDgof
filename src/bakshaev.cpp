#include <Rcpp.h>
#include <cmath> // For std::sqrt, std::max
using namespace Rcpp;

//' Find evaluation points
//'
//' @param rnull, a function that generate new data.
//' @param p, a vector of parameters for rnull.
//' @param m, size of matrix.
//' @return a matrix
//' 
// [[Rcpp::export]]
NumericMatrix gen_eval(Function rnull,
                       NumericVector p,
                       int m) {
  int i, j, k, n, d;
  NumericMatrix simdta;
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  if(resrnull.size()==0) simdta = rnull();
  else simdta = rnull(p);
  n=simdta.nrow(), d=simdta.ncol();
  NumericMatrix Eval(m, d);
  k=-1;
  do {
    ++k;
    if(k>0) {
      if(resrnull.size()==0) simdta = rnull();
      else simdta = rnull(p);
    }  
    for(i=0;i<n;++i) {
      if(i+k*n<m) {
        for(j=0;j<d;++j)
          Eval(i+k*n,j)=simdta(i,j);
      }
    }
  } while ((k+1)*n<m); 
  return Eval;
}  
   
//' Find gaussian kernel pdf
//' 
//' @param Eval a matrix.
//' @param S a matrix
//' @param h bandwith, a double
//' @return a matrix
//' 
// [[Rcpp::export]]
NumericMatrix gauss_kernel_matrix(const NumericMatrix& Eval,
                                  const NumericMatrix& S,
                                  double h) {
  int m = Eval.nrow();
  int d = Eval.ncol();
  int n = S.nrow();
  
  NumericMatrix K(m, n);
  double norm_const = std::pow(2.0 * M_PI, -0.5 * d) * std::pow(h, -d);
  double inv_2h2 = 1.0 / (2.0 * h * h);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double dist2 = 0.0;
      for (int k = 0; k < d; k++) {
        double diff = Eval(i, k) - S(j, k);
        dist2 += diff * diff;
      }
      K(i, j) = norm_const * std::exp(-dist2 * inv_2h2);
    }
  }
  return K;
}

//' Estimate E and Var/n at Eval for given h, using MC from rnull
//' 
//' @param rnull generate data under the null hypothesis.
//' @param p values for rnull
//' @param Eval matrix of evaluations
//' @param h bandwith, a double
//' @param nsim_mc number of simulation runs
//' @param n sample size
//' @return a matrix
//'
// [[Rcpp::export]]
List estimateEV(Function rnull,
                            NumericVector p,
                            const NumericMatrix& Eval,
                            double h,
                            int nsim_mc,
                            int n) {
  int m = Eval.nrow(), d = Eval.ncol();
  NumericMatrix Y(nsim_mc, d), simdta(n, d);
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  Y=gen_eval(rnull, p, nsim_mc);
  NumericMatrix Kmc = gauss_kernel_matrix(Eval, Y, h);
  
  NumericVector Ekh(m);
  for (int i = 0; i < m; i++) {
    double sum = 0.0;
    for (int t = 0; t < nsim_mc; t++) sum += Kmc(i, t);
    Ekh[i] = sum / nsim_mc;
  }
  
  NumericVector Var_fhat(m);
  for (int i = 0; i < m; i++) {
    double acc = 0.0;
    for (int t = 0; t < nsim_mc; t++) {
      double diff = Kmc(i, t) - Ekh[i];
      acc += diff * diff;
    }
    double var_pop = acc / nsim_mc;
    Var_fhat[i] = var_pop / n;
  }
  
// adaptive floor
  double median_var = 0.0;
  std::vector<double> positives;
  for (int i = 0; i < m; i++) if (Var_fhat[i] > 0 && R_finite(Var_fhat[i])) positives.push_back(Var_fhat[i]);
  if (!positives.empty()) {
    std::nth_element(positives.begin(), positives.begin() + positives.size()/2, positives.end());
    median_var = positives[positives.size()/2];
  }
  double floor_var = (median_var <= 0.0) ? 1e-12 : median_var * 1e-6;
  for (int i = 0; i < m; i++) {
    if (!R_finite(Var_fhat[i]) || Var_fhat[i] < floor_var) Var_fhat[i] = floor_var;
  }
  
  return List::create(Named("Ekh") = Ekh,
                      Named("Var_fhat") = Var_fhat,
                      Named("floor_var") = floor_var);
}

//' Compute M statistic for one dtaset
//' 
//' @param dta data matrix.
//' @param Eval matrix of evaluations
//' @param hs bandwiths
//' @param rnull generate new data
//' @param p values for parametric bootstrap
//' @param nsim_mc number of simulation runs
//' @return a double
//'
// [[Rcpp::export]]
double compute_M_for_dtaset(const NumericMatrix& dta,
                             const NumericMatrix& Eval,
                             const NumericVector& hs,
                             Function rnull,
                             NumericVector p, 
                             int nsim_mc) {
  int m = Eval.nrow();
  double M_val = -INFINITY;
  
  for (int j = 0; j < hs.size(); j++) {
    double h = hs[j];
    
    List EV = estimateEV(rnull, p, Eval, h, nsim_mc, dta.nrow());
    NumericVector Ekh = EV["Ekh"];
    NumericVector Var_fhat = EV["Var_fhat"];
    
    NumericMatrix Kobs = gauss_kernel_matrix(Eval, dta, h);
    
    NumericVector fhat(m);
    for (int i = 0; i < m; i++) {
      double sum = 0.0;
      for (int t = 0; t < dta.nrow(); t++) sum += Kobs(i, t);
      fhat[i] = sum / dta.nrow();
    }
    
    NumericVector xi(m);
    for (int i = 0; i < m; i++) {
      double denom = std::sqrt(Var_fhat[i]);
      if (!R_finite(denom) || denom <= 0.0) denom = 1e-8;
      xi[i] = (fhat[i] - Ekh[i]) / denom;
      if (!R_finite(xi[i])) xi[i] = 0.0;
    }
    
    double maxabs = 0.0;
    for (int i = 0; i < m; i++) {
      double val = std::abs(xi[i]);
      if (val > maxabs) maxabs = val;
    }
    if (maxabs > M_val) M_val = maxabs;
  }
  return M_val;
}
//' Run Bakshaev and Rudzkis Test
//' 
//' @param dta data matrix.
//' @param rnull generate new data.
//' @param p , parameters for parametric bootstrap.
//' @param m_eval =100, number of evaluation points of kde.
//' @param nsim =200, number of simulation runs.
//' @param nsim_mc =1000, number of simulation runs.
//' @return a list
//'
// [[Rcpp::export]]
List bakshaev_rudzkis(NumericMatrix dta,
                      Function rnull,
                      NumericVector p,
                      int m_eval = 100,
                      int nsim = 200,
                      int nsim_mc = 1000) {
  int n = dta.nrow(), d = dta.ncol();
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  NumericMatrix Eval(m_eval,d), simdta(n, d);
  NumericVector hs(3);
  hs(0)=0.3;hs(1)=0.5;hs(2)=0.8;
  Eval=gen_eval(rnull, p, m_eval);
  double M_obs = compute_M_for_dtaset(dta, Eval, hs, rnull, p, nsim_mc);
  NumericVector null_stats(nsim);
  for (int s = 0; s < nsim; s++) {
    if(resrnull.size()==0) simdta = rnull();
    else simdta = rnull(p);
    null_stats[s] = compute_M_for_dtaset(simdta, Eval, hs, rnull, p,nsim_mc);
  }
  
  double pval = 0.0;
  for (int s = 0; s < nsim; s++) if (null_stats[s] >= M_obs) pval += 1.0;
  pval /= nsim;
  
  return List::create(Named("statistic") = M_obs,
                      Named("p.value") = pval);
}
