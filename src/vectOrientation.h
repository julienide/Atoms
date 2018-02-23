#ifndef VECT_ORIENTATION_H
#define VECT_ORIENTATION_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

Rcpp::NumericVector vectOrientation(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1) ;

Rcpp::NumericVector vectOrientationAxis(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::NumericVector V2) ;

Rcpp::NumericVector vectOrientationPairwise(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::DataFrame& V2) ;

Rcpp::NumericVector vectOrientationNotPairwise(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::DataFrame& V2) ;


// 
// // Angles between a set of vectors and the x,y,z-axis
// Rcpp::S4 vectOrientationXYZ(
//     const Rcpp::S4& x,
//     const Rcpp::DataFrame& V1,
//     const Rcpp::IntegerVector& V2) ;
// 
// // Angles between two sets of plan normals
// Rcpp::S4 planOrientation(
//     const Rcpp::S4& x,
//     const Rcpp::DataFrame& V1,
//     const Rcpp::DataFrame& V2) ;
// 
// // Angles between a set of plan normals and the x,y,z-axis
// Rcpp::S4 planOrientationXYZ(
//     const Rcpp::S4& x,
//     const Rcpp::DataFrame& V1,
//     const Rcpp::IntegerVector& V2) ;

#endif // VECT_ORIENTATION_H
