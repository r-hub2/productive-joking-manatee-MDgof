#include <Rcpp.h>
using namespace Rcpp;
//' This function performs a Rosenblatt transform
//' @param  x data set (a matrix)
//' @param  cdf distribution function
//' @param  p (possible) parameters for cdf
//' @param  Range matrix with range of data
//' @keywords internal
//' @return A matrix of transformed data
// [[Rcpp::export]]
NumericMatrix rosenblattC(NumericMatrix x, Function cdf, 
                          NumericVector p,
                          NumericMatrix Range) {
   Rcpp::Environment base("package:base");
   Rcpp::Function formals_r = base["formals"];
   Rcpp::List rcdf = formals_r(Rcpp::_["fun"]=cdf);
   int i, n=x.nrow();
   double h=(max(x(_, 0))-min(x(_, 0)))/1000000.0;
   double Ry=Range(1, 1);
   if (!(R_finite(Ry))) Ry=2*max(x(_, 1));
   
   NumericMatrix x1(n,2),x2(n,2),x3(n,2),x4(n,2), x5(n,2), y(n,2);
   NumericVector u(n), num1(n), num2(n), den1(n), den2(n);
   for(i=0;i<n;++i) {
     x1(i, 0)=x(i,0);
     x1(i, 1)=Ry;
     x2(i, 0)=x(i,0);
     if(x(i,0)<Range(1, 0)-h) x2(i, 0)=x(i,0)+h;
     x2(i, 1)=Ry;
     x3(i, 0)=x(i,0);
     if(x(i,0)>Range(0, 0)+h) x3(i, 0)=x(i,0)-h;
     x3(i, 1)=Ry;
     x4(i, 0)=x(i,0);
     if(x(i,0)<Range(1, 0)-h) x4(i, 0)=x(i,0)+h;
     x4(i, 1)=x(i,1);
     x5(i, 0)=x(i,0);
     if(x1(i,0)<Range(0, 0)+h) x5(i, 0)=x(i,0)-h;
     x5(i, 1)=x(i,1);
   }
   if(rcdf.size()==1) {
      u=cdf(x1);
      num1=cdf(x4);
      num2=cdf(x5);
      den1=cdf(x2);
      den2=cdf(x3);
   }
   else {
     u=cdf(x1, p);
     num1=cdf(x4, p);
     num2=cdf(x5, p);
     den1=cdf(x2, p);
     den2=cdf(x3, p);
   }
   y(_,0)=u;
   for(i=0;i<n;++i) {
      if(std::abs(den1[i]-den2[i])<1e-12) {
         y(i,1)=runif(1, 0, 1e-5)[0];
      }  
      else y(i,1)=(num1[i]-num2[i])/(den1[i]-den2[i]);
   }    
   return y;
 }
