#include "LAMMPSParser.h"

//' @export
//' @rdname Atoms
//' @noRd
// [[Rcpp::export(name = ".LAMMPSParser")]]
Rcpp::S4 LAMMPSParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector first = Rcpp::IntegerVector::create(1),
    Rcpp::IntegerVector last = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1)){
  
  LAMMPSReader obj(files, selection, first, last, stride) ;
  obj.readAllFiles() ;
  return obj.asS4() ;
}