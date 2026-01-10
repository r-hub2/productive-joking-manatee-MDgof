#include <Rcpp.h>
using namespace Rcpp;

//' Bins continuous data
//' 
//' @param x A numeric matrix with two columns.
//' @param nbins number of bins.
//' @param Range range of variables
//' @param ChangeVals =FALSE, should values of discrete rv's be adjusted to midpoints?
//' @return A numeric matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix discretize(
        Rcpp::NumericMatrix x, 
        Rcpp::NumericMatrix Range,
        Rcpp::NumericVector nbins,
        bool ChangeVals=false) {
  Range.attr("dim") = Dimension(2, 2);
  int n=x.nrow(), i, j, k;
  NumericMatrix grid(nbins[0]*nbins[1], 3);
  k=0;
  double mx, my, Mx, My;
  LogicalVector tmp=is_infinite(Range);
  if(tmp(0)) mx=min(x(_,0));
  else mx=Range(0,0);
  if(tmp(1)) Mx=max(x(_,0));
  else Mx=Range(1,0);
  if(tmp(2)) my=min(x(_,1));
  else my=Range(0,1);
  if(tmp(3)) My=max(x(_,1));  
  else My=Range(1,1);
  double stepx=(Mx-mx)/nbins[0],stepy=(My-my)/nbins[1];
  for(j=0;j<nbins[1];++j) {
    for(i=0;i<nbins[0];++i) {
       grid(k, 0) = mx+(i+1)*stepx;
       grid(k, 1) = my+(j+1)*stepy;
       ++k;
    }
  }
  for(k=0;k<n;++k) {
    i=0;
    while( ((x(k,0)>grid(i,0)) || 
           (x(k,1)>grid(i,1)) )&& 
           (i<nbins[0]*nbins[1]-1) ) ++i;
    grid(i, 2)=grid(i, 2)+1.0;
  }
  if(ChangeVals) {
     for(i=0;i<grid.nrow();++i) {
       grid(i,0)=grid(i,0)-stepx/2.0; 
       grid(i,1)=grid(i,1)-stepy/2.0;
     }   
  }
  
  return grid;
} 
