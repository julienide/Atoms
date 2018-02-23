#include "PDBParser.h"

//' @export
//' @rdname Atoms
//' @noRd
// [[Rcpp::export(name = ".PDBParser")]]
Rcpp::S4 PDBParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector first = Rcpp::IntegerVector::create(1),
    Rcpp::IntegerVector last = Rcpp::IntegerVector::create(-1),
    Rcpp::IntegerVector stride = Rcpp::IntegerVector::create(1)){
  
  PDBReader obj(files, selection, first, last, stride) ;
  obj.readAllFiles() ;
  return obj.asS4() ;
}