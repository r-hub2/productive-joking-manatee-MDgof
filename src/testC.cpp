#include <Rcpp.h>
using namespace Rcpp;
#include "calcTS.h"
//' run gof tests for continuous data
//' 
//' @param dta A numeric matrix of data
//' @param rnull R function (generate data under null hypothesis)
//' @param TS function that calculates test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B (=5000) Number of simulation runs 
//' @keywords internal
//' @return A matrix of numbers (test statistics and p values)
// [[Rcpp::export]]
Rcpp::NumericMatrix testC(
        Rcpp::NumericMatrix dta, 
        Rcpp::Function rnull, 
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        int B=5000) {

  int i,j;
  NumericVector TS_data=calcTS(dta, TS, typeTS, TSextra);
  
  int const nummethods=TS_data.size();
  Rcpp::CharacterVector methods=TS_data.names();
  NumericVector dtasim, TS_sim(nummethods), pvals(nummethods), p;
  NumericMatrix out(2, nummethods);
  colnames(out) = methods;
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  Rcpp::Function phat=TSextra["phat"];
  p=phat(dta);
  for(i=0;i<B;++i) {
    if(resrnull.size()==0) dta=rnull();
    else dta=rnull(p);
    NumericVector TS_sim=calcTS(dta, TS, typeTS, TSextra);
    for(j=0;j<nummethods;++j) {
      if(TS_data(j)<TS_sim(j)) pvals(j)=pvals(j)+1;
    }
  }
  for(i=0;i<nummethods;++i) {
      out(0, i)=TS_data(i);
      out(1, i)=pvals(i)/B;
  }
  return out;
}
