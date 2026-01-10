#ifndef BAKSHAEV_RUDZKIS_H
#define BAKSHAEV_RUDZKIS_H

#include <Rcpp.h> 
Rcpp::List bakshaev_rudzkis(Rcpp::NumericMatrix dta,
                       Rcpp::Function rnull,
                       Rcpp::NumericVector p,
                       int m_eval = 100,
                       int nsim = 200,
                       int nsim_mc = 1000);
#endif

#ifndef GAUSS_KERNEL_MATRIX_H
#define GAUSS_KERNEL_MATRIX_H

#include <Rcpp.h>
Rcpp::NumericMatrix gauss_kernel_matrix(
              const Rcpp::NumericMatrix& Eval,
              const Rcpp::NumericMatrix& S,
              double h);
#endif

#ifndef GEN_EVAL_H
#define GEN_EVAL_H

#include <Rcpp.h>
Rcpp::NumericMatrix gen_eval(
              Rcpp::Function rnull,
              Rcpp::NumericVector p,
              int m);
#endif  

#ifndef ESTIMATE_E_AND_VARFHAT_H
#define ESTIMATE_E_AND_VARFHAT_H

#include <Rcpp.h>
Rcpp::List estimate_E_and_Varfhat(
              Rcpp::Function rnull,
              Rcpp::NumericVector p,
              const Rcpp::NumericMatrix& Eval,
              Rcpp::NumericMatrix S,
              double h,
              int nsim_mc,
              int n);
#endif

#ifndef COMPUTE_M_FOR_DTASET_H
#define COMPUTE_M_FOR_DTASET_H

#include <Rcpp.h>
double compute_M_for_dtaset(const Rcpp::NumericMatrix& dta,
                             const Rcpp::NumericMatrix& Eval,
                             const Rcpp::NumericVector& hs,
                             Rcpp::Function rnull,
                             Rcpp::NumericVector p, 
                             int nsim_mc);
#endif  

                             
