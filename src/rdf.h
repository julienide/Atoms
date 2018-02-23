#ifndef RDF_H
#define RDF_H

// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>
#include "applyPBC.h"

Rcpp::DataFrame rdfAtoms(
    const Rcpp::S4& x,
    const std::vector< int >& sel1,
    const std::vector< int >& sel2,
    const std::vector< int >& resnumb,
    const double& cutoff,
    const double& interval,
    const std::vector< bool >& type) ;

#endif // RDF_H
