#include <Rcpp.h>
#include "bakshaev.h"
#include "p2dC.h"
using namespace Rcpp;

//' Find test statistics for discrete data
//' 
//' @param x A numeric matrix.
//' @param pnull cdf.
//' @param param parameters for pnull  in case of parameter estimation.
//' @param TSextra list with additional info
//' @return A numeric vector with test statistics
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector TS_disc(
        Rcpp::NumericMatrix x, 
        Rcpp::Function pnull,
        Rcpp::NumericVector param,
        Rcpp::List TSextra) {
  Rcpp::CharacterVector methods=Rcpp::CharacterVector::create("KS", "K", "CvM", "AD", "TV", "KL", "H", "P");
  int const nummethods=methods.size();
  int n=x.nrow(), i;
  NumericVector TS(nummethods);
  double tmp;
  TS.names() =  methods;
  /* some setup */
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List respnull = formals_r(Rcpp::_["fun"]=pnull);
  NumericVector Fx(n);
  NumericMatrix valsx=x(_, Range(0, 1));
  if(respnull.size()==1) Fx = pnull(valsx);
  else Fx = pnull(valsx, param);
 /*density*/
  NumericVector dnull=p2dC(x, pnull, param, Fx)(_,4);
  double M=0.0;
  for(i=0;i<n;++i) M=M+x(i,2);
  NumericVector ecdf(n);
  for(i=0;i<n;++i) {
    for(int j=0;j<n;++j) {
      if((x(j,0)<=x(i,0)) && (x(j,1)<=x(i,1)))
         ecdf(i)=ecdf(i)+x(j,2);
    }  
    ecdf(i)=ecdf(i)/M;
  }
  /*  Kolmogorov-Smirnov and Kuiper*/
  TS(0)=0;
  double mx=0.0, Mx=0.0;
  for(i=0;i<n;++i) {
    tmp = Fx[i]-ecdf(i);
    if(tmp<mx) mx=tmp;
    if(tmp>Mx) Mx=tmp;
  }
  if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
  else TS(0)=std::abs(Mx);
  TS(1) = std::abs(mx)+std::abs(Mx);
    
  /* Cramer-vonMises and Anderson-Darling */
  
  TS(2)=0;
  TS(3)=0;
  for(i=0;i<n;++i) {
    tmp = Fx[i]-ecdf(i);
    TS(2) = TS(2) + tmp*tmp*dnull(i);
    if((Fx[i]>0) && (Fx[i]<1))
      TS(3) = TS(3) + tmp*tmp/Fx[i]/(1-Fx[i])*dnull(i);
  }
  
  /* Total Variation, Kullback-Leibler, Hellinger, Pearson*/
  
  for(i=0;i<n;++i) {
    double on=x(i, 2)/M;
    TS(4)=TS(4)+std::abs(on-dnull(i));
    if(on>0) TS(5)=TS(5)+on*log(on/dnull(i));
    tmp=std::sqrt(on)-std::sqrt(dnull(i));
    TS(6)=TS(6)+tmp*tmp;
    tmp=dnull(i);
    if(tmp>0) TS(7)=TS(7)+M*(on-tmp)*(on-tmp)/tmp;
  }
  return TS;
} 
