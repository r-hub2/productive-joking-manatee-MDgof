#include <Rcpp.h>
#include "calcTS.h"
using namespace Rcpp;

//' find power of gof tests for continuous data
//' 
//' @param rnull R function (generate data under null hypothesis)
//' @param ralt  R function to generate data under alternative
//' @param param_alt parameters of ralt
//' @param TS function to calculate test statistics
//' @param typeTS integer indicating type of test statistic
//' @param TSextra list to pass to TS
//' @param B  =1000 Number of simulation runs
//' @keywords internal
//' @return A matrix of powers
// [[Rcpp::export]]
Rcpp::List powerC(
        Rcpp::Function rnull,
        Rcpp::Function ralt, 
        Rcpp::NumericVector param_alt,
        Rcpp::Function TS,
        int typeTS,
        Rcpp::List TSextra,
        const int B=1000) {

  int  i, l, m, cn, np=param_alt.size();
  NumericMatrix dta=ralt(param_alt[0]);
  NumericVector TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List resrnull = formals_r(Rcpp::_["fun"]=rnull);
  Rcpp::List resralt = formals_r(Rcpp::_["fun"]=ralt);
  Rcpp::Function phat=TSextra["phat"];
  Rcpp::Function knn=TSextra["knn"];
  NumericVector p=phat(dta), param_alt_vec(B*np);
  NumericMatrix realdata(B, nummethods), simdata(B*np, nummethods);
  Rcpp::CharacterVector methods=TS_data.names();
  cn=-1;
  for(l=0;l<B;++l) {
     if(resrnull.size()==0) dta=rnull();
     else dta=rnull(p);
     if(!TSextra["NoDensity"]) TSextra["nn"]=knn(dta);
     TS_data = calcTS(dta, TS, typeTS, TSextra); 
     for(i=0;i<nummethods;++i) realdata(l,i)=TS_data(i);
     for(m=0;m<np;++m) {
         ++cn;
         param_alt_vec[cn]=param_alt[m];
         if(resralt.size()==1) dta=ralt(param_alt[m]);
         else dta=ralt(param_alt[m], p);
         if(!TSextra["NoDensity"]) TSextra["nn"]=knn(dta);
         TS_data = calcTS(dta, TS, typeTS, TSextra);
         for(i=0;i<nummethods;++i) simdata(cn, i)=TS_data(i);
     }   
  }
  return List::create(Named("Data")=realdata, 
                      Named("Sim")=simdata,
                      Named("param_alt_vec")=param_alt_vec);
}
