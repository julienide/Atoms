#ifndef PDB_PARSER_H
#define PDB_PARSER_H

#include <Rcpp.h>
#include "PDBReader.h"

Rcpp::S4 PDBParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection,
    Rcpp::IntegerVector first,
    Rcpp::IntegerVector last,
    Rcpp::IntegerVector stride) ;

#endif // PDB_PARSER_H
