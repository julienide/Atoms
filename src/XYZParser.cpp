#include "XYZParser.h"

//' @export
//' @rdname Atoms
//' @noRd
// [[Rcpp::export(name = ".XYZParser")]]
Rcpp::S4 XYZParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector first = Rcpp::IntegerVector::create(1),
    Rcpp::IntegerVector last = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1)){
  
  XYZReader obj(files, selection, first, last, stride) ;
  obj.readAllFiles() ;
  return obj.asS4() ;
}