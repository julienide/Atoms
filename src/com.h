#ifndef COM_H
#define COM_H

#include <Rcpp.h>

Rcpp::S4 comAtoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable<Rcpp::NumericVector>& mass,
    const Rcpp::Nullable<Rcpp::IntegerVector>& factor) ;

#endif // COM_H
