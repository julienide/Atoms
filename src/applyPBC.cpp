#include "applyPBC.h"

void applyPBC(
    std::vector< double >& d,
    const Rcpp::NumericMatrix& cell,
    const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc){
  std::vector< double > dabc(3) ;
  dabc[0] = d[0]*cellInv(0, 0) + d[1]*cellInv(0, 1) + d[2]*cellInv(0, 2) ;
  dabc[1] = d[0]*cellInv(1, 0) + d[1]*cellInv(1, 1) + d[2]*cellInv(1, 2) ;
  dabc[2] = d[0]*cellInv(2, 0) + d[1]*cellInv(2, 1) + d[2]*cellInv(2, 2) ;
  for(int i = 0; i < 3; i++){
    if(pbc[i]){
      dabc[i] = dabc[i] - floor(dabc[i] + 0.5) ;
    }
  }
  d[0] = dabc[0]*cell(0, 0) + dabc[1]*cell(0, 1) + dabc[2]*cell(0, 2) ;
  d[1] = dabc[0]*cell(1, 0) + dabc[1]*cell(1, 1) + dabc[2]*cell(1, 2) ;
  d[2] = dabc[0]*cell(2, 0) + dabc[1]*cell(2, 1) + dabc[2]*cell(2, 2) ;
}

// This version should be mush faster
void applyPBC2(
    double& dx, double& dy, double& dz,
    const Rcpp::NumericMatrix& cell,
    const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc){
  double da = dx*cellInv(0, 0) + dy*cellInv(0, 1) + dz*cellInv(0, 2) ;
  double db = dx*cellInv(1, 0) + dy*cellInv(1, 1) + dz*cellInv(1, 2) ;
  double dc = dx*cellInv(2, 0) + dy*cellInv(2, 1) + dz*cellInv(2, 2) ;
  if(pbc[0]) {
    da = da - floor(da + 0.5) ;
  }
  if(pbc[1]) {
    db = db - floor(db + 0.5) ;
  }
  if(pbc[2]) {
    dc = dc - floor(dc + 0.5) ;
  }
  dx = da*cell(0, 0) + db*cell(0, 1) + dc*cell(0, 2) ;
  dy = da*cell(1, 0) + db*cell(1, 1) + dc*cell(1, 2) ;
  dz = da*cell(2, 0) + db*cell(2, 1) + dc*cell(2, 2) ;
}