#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics for  data
//' @param  dta data set (a matrix)
//' @param  TS routine
//' @param  typeTS format of TS
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of values of test statistic(s)
// [[Rcpp::export]]
NumericVector calcTS(
   Rcpp::NumericMatrix dta, 
   Rcpp::Function TS,
   int typeTS,
   Rcpp::List TSextra) {
   NumericVector TS_data;
   Function phat=TSextra["phat"];
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   NumericVector p=phat(dta);
   if(typeTS==1) TS_data=TS(dta, TSextra["pnull"], p);
   if(typeTS==2) TS_data=TS(dta, TSextra["pnull"], p, TSextra);
   return  TS_data;
 }
 
 
