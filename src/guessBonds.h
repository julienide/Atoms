#ifndef GUESSBONDS_H
#define GUESSBONDS_H

#include <Rcpp.h>

Rcpp::DataFrame guessBondsAtoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable< Rcpp::NumericVector > radius,
    const Rcpp::NumericVector& safety,
    const Rcpp::LogicalVector& usePBC) ;

#endif // GUESSBONDS_H
