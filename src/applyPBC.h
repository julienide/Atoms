#include <Rcpp.h>

void applyPBC(
    std::vector< double >& d,
    const Rcpp::NumericMatrix& cell,
    const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc) ;

void applyPBC2(
    double& dx, double& dy, double& dz,
    const Rcpp::NumericMatrix& cell,
    const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc) ;