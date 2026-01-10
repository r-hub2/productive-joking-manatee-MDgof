#ifndef CALCTS_H
#define CALCTS_H

#include <Rcpp.h>
Rcpp::NumericVector calcTS(
       Rcpp::NumericMatrix dta, 
       Rcpp::Function TS, 
       int typeTS,
       Rcpp::List TSextra
        );

#endif
