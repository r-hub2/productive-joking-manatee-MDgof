#include <Rcpp.h>
#include "mdecdf.h"
#include "bakshaev.h"
#include "mmd.h"
#include "ksd.h"
#include "rosenblattC.h"
#include "FF.h"
using namespace Rcpp;

//' Find test statistics for continuous data
//' 
//' @param x A numeric matrix.
//' @param pnull cdf.
//' @param param parameters for pnull  in case of parameter estimation.
//' @param TSextra list with additional info
//' @return A numeric vector with test statistics
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector TS_cont(
        Rcpp::NumericMatrix x, 
        Rcpp::Function pnull,
        Rcpp::NumericVector param,
        Rcpp::List TSextra) {
  Rcpp::CharacterVector methods=Rcpp::CharacterVector::create("qKS", "qK", "qCvM", "qAD", "BB", "BR", "MMD", "KSD", "FF", "RK");
  int const nummethods=methods.size();
  int n=x.nrow(), Dim=x.ncol(), i, j;
  NumericVector TS(nummethods);
  double tmp;
  TS.names() =  methods;

  /* some setup */
  Rcpp::Environment base("package:base");
  Rcpp::Function formals_r = base["formals"];
  Rcpp::List respnull = formals_r(Rcpp::_["fun"]=pnull);
  NumericVector Fx(n);
  if(respnull.size()==1) Fx = pnull(x);
  else Fx = pnull(x, param);

  /*  empirical distribution function */
  
  NumericVector ecdf=mdecdf(x, x);

  /*  Kolmogorov-Smirnov and Kuiper*/
  TS(0)=0;
  double mx=0.0, Mx=0.0;
  for(i=0;i<n;++i) {
    tmp = Fx[i]-ecdf[i];
    if(tmp<mx) mx=tmp;
    if(tmp>Mx) Mx=tmp;
    tmp = Fx[i]-ecdf[i]+1.0/n;
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
    tmp = Fx[i]-ecdf[i];
    TS(2) = TS(2) + tmp*tmp;
    if((Fx[i]>0) && (Fx[i]<1))
      TS(3) = TS(3) + tmp*tmp/Fx[i]/(1-Fx[i]);
  }
  
  /* Bickel-Breiman, if density was supplied*/
  
  if(!TSextra["NoDensity"]) {
    NumericVector fx(n);
    Function dnull=TSextra["dnull"];
    if(respnull.size()==1) fx = dnull(x);
    else fx = dnull(x, param);
    
    Rcpp::List nn = TSextra["nn"];
    NumericMatrix nn1=nn["nn.dist"];
    NumericVector R=nn1(_,0);
    NumericVector U(n);
    double K=3.141593;
    if(Dim>2) 
      for(j=3;j<=Dim;++j) K=K*1.772454;
    tmp=::tgamma(Dim/2.0+1.0);
    K=K/tmp;
    TS[4]=0;
    for(i=0;i<n;++i) {
      tmp=R[i]*R[i];
      if(Dim>2) 
        for(j=3;j<=Dim;++j) tmp=tmp*R[i];
      U[i]=exp(-n*fx[i]*K*tmp);
    }
    std::sort(U.begin(), U.end());
    for(i=0;i<n;++i) 
       TS[4]=TS[4]+(U[i]-(i+1.0)/n)*(U[i]-(i+1.0)/n);
    TS[4]=TS[4]/n;
    
  }
 
  /* Bakshaev-Rudzkis*/ 
  
  Function rnull=TSextra["rnull"];
  NumericVector hs(3);
  hs(0)=0.3;hs(1)=0.5;hs(2)=0.8;
  NumericMatrix Eval=TSextra["Eval"];
  TS[5]=compute_M_for_dtaset(x, Eval, hs, rnull, param, 1000);
    
  /* Maximum Mean Discrepancy */
  
  TS[6]=mmd(x, TSextra["rnull"], param);
  if(!TSextra["NoDensity"])
     TS[7]=ksd(x, TSextra["scf"], param);
  /* Fasano-Franceschini and Ripley's K tests */
  
  if(Dim==2) {
     NumericMatrix transx=rosenblattC(x, pnull, param, TSextra["Range"]);
     TS[8]=FF(transx);
     Function rk=TSextra["ripleyK"];
     NumericVector tmprk=rk(transx);
     TS[9]=tmprk[0];
  }  
  return TS;
} 
