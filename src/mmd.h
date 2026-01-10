#ifndef MMD_H
#define MMD_H

#include <Rcpp.h>
double mmd(Rcpp::NumericMatrix x,
           Rcpp::Function rnull,
           Rcpp::NumericVector p);
#endif
