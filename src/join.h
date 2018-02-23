#ifndef JOIN_H
#define JOIN_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

void joinAtoms(
    Rcpp::S4& x,
    const Rcpp::LogicalVector& multi) ;

#endif // JOIN_H