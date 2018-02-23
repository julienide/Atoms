#ifndef FREE_VOLUME_H
#define FREE_VOLUME_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

Rcpp::NumericVector freeVolumeAtoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable< Rcpp::NumericVector > radius) ;

Rcpp::NumericVector freeVolume2Atoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable< Rcpp::NumericVector > radius) ;
  
#endif // FREE_VOLUME_H
