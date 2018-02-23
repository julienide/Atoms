#ifndef PLAN_NORMALS_H
#define PLAN_NORMALS_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

Rcpp::S4 planNormalsAtoms(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3) ;

#endif // PLAN_NORMALS_H