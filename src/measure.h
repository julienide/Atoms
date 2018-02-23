#ifndef MEASURE_H
#define MEASURE_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

Rcpp::NumericMatrix measureBonds(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2) ;

Rcpp::NumericMatrix measureAngles(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3) ;

Rcpp::NumericMatrix measureDihedrals(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3,
    const Rcpp::IntegerVector& atm4) ;

#endif // MEASURE_H