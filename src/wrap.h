#ifndef WRAP_H
#define WRAP_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>

void wrapAtoms(
    Rcpp::S4& x,
    const Rcpp::IntegerVector& resnumb,
    const Rcpp::LogicalVector& multi) ;

#endif // WRAP_H