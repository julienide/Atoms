#ifndef MOL2_PARSER_H
#define MOL2_PARSER_H

#include <Rcpp.h>
#include "MOL2Reader.h"

Rcpp::S4 MOL2Parser(
  Rcpp::CharacterVector files,
  Rcpp::IntegerVector selection,
  Rcpp::IntegerVector first,
  Rcpp::IntegerVector last,
  Rcpp::IntegerVector stride) ;

#endif // MOL2_PARSER_H
