#include "MOL2Parser.h"

//' @export
//' @rdname Atoms
//' @noRd
// [[Rcpp::export(name = ".MOL2Parser")]]
Rcpp::S4 MOL2Parser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector first = Rcpp::IntegerVector::create(1),
    Rcpp::IntegerVector last = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1)){
  
  MOL2Reader obj(files, selection, first, last, stride) ;
  obj.readAllFiles() ;
  return obj.asS4() ;
}