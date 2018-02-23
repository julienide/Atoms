#ifndef RMSD_H
#define RMSD_H

#include <Rcpp.h>
#include "applyPBC.h"

Rcpp::NumericVector rmsdAtoms(
    const Rcpp::S4& x,
    const Rcpp::S4& y,
    const Rcpp::LogicalVector& usePBC) ;

#endif // RMSD_H
