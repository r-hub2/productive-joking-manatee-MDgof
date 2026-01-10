#ifndef KSD_H
#define KSD_H

#include <Rcpp.h>
double ksd(
       Rcpp::NumericMatrix x, 
       Rcpp::Function scf, 
       Rcpp::NumericVector p
        );

#endif
