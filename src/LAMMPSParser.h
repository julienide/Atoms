#ifndef LAMMPS_PARSER_H
#define LAMMPS_PARSER_H

#include <Rcpp.h>
#include "LAMMPSReader.h"

Rcpp::S4 LAMMPSParser(
    Rcpp::CharacterVector files,
    Rcpp::IntegerVector selection,
    Rcpp::IntegerVector first,
    Rcpp::IntegerVector last,
    Rcpp::IntegerVector stride) ;

#endif // LAMMPS_PARSER_H
