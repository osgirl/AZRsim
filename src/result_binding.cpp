#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
// Function to replace rbind function of R
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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
