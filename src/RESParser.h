#ifndef RES_PARSER_H
#define RES_PARSER_H

#include <Rcpp.h>
#include "RESReader.h"

Rcpp::S4 RESParser(
  Rcpp::CharacterVector files,
  Rcpp::IntegerVector selection,
  Rcpp::IntegerVector first,
  Rcpp::IntegerVector last,
  Rcpp::IntegerVector stride) ;

#endif // RES_PARSER_H
