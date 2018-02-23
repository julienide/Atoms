#ifndef XYZ_PARSER_H
#define XYZ_PARSER_H

#include <Rcpp.h>
#include "XYZReader.h"

Rcpp::S4 XYZParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection,
    Rcpp::IntegerVector first,
    Rcpp::IntegerVector last,
    Rcpp::IntegerVector stride) ;

#endif // XYZ_PARSER_H
