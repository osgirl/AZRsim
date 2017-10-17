#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix result_binding(NumericMatrix simres_all, NumericMatrix simres_piece)
{
  int result_nrow = simres_all.nrow()+simres_piece.nrow();
  int result_ncol = simres_all.ncol();
  NumericMatrix result(result_nrow,result_ncol);
  int k=0;
  for(int i=0;i<simres_all.nrow();i++)
  {
    result(k++,_)=simres_all(i,_);
  }
  for(int i=0;i<simres_piece.nrow();i++)
  {
    result(k++,_)=simres_piece(i,_);
  }
  return result;

}
