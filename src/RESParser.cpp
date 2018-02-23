#include "RESParser.h"

//' @export
//' @rdname Atoms
//' @noRd
// [[Rcpp::export(name = ".RESParser")]]
Rcpp::S4 RESParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector first = Rcpp::IntegerVector::create(1),
    Rcpp::IntegerVector last = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1)){
  
  RESReader obj(files, selection, first, last, stride) ;
  obj.readAllFiles() ;
  return obj.asS4() ;
}