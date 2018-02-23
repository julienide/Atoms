#include "rmsd.h"

//' @rdname rmsd
//' @noRd
//' @export
// [[Rcpp::export(name = ".rmsdAtoms")]]
Rcpp::NumericVector rmsdAtoms(
    const Rcpp::S4& x,
    const Rcpp::S4& y,
    const Rcpp::LogicalVector& usePBC)
{
  Rcpp::NumericMatrix xPx = x.attr("x") ;
  Rcpp::NumericMatrix xPy = x.attr("y") ;
  Rcpp::NumericMatrix xPz = x.attr("z") ;
  Rcpp::LogicalVector xPBC = x.attr("pbc") ;
  Rcpp::NumericMatrix xA = x.attr("a") ;
  Rcpp::NumericMatrix xB = x.attr("b") ;
  Rcpp::NumericMatrix xC = x.attr("c") ;
  int nframe = xPx.ncol() ;
  int natom = xPx.nrow() ;
  
  Rcpp::NumericMatrix yPx = y.attr("x") ;
  Rcpp::NumericMatrix yPy = y.attr("y") ;
  Rcpp::NumericMatrix yPz = y.attr("z") ;

  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  Rcpp::NumericVector rmsd(nframe) ;
  int n = 0 ;
  for(int frame = 0; frame < nframe; frame++){
    Rcpp::NumericMatrix cell(3, 3) ;
    Rcpp::NumericMatrix cellInv(3, 3) ;
    if(usePBC[0] && (xPBC[0] || xPBC[1] || xPBC[2])){
      cell(0, 0) = xA(0, frame) ; cell(0, 1) = xB(0, frame) ; cell(0, 2) = xC(0, frame) ;
      cell(1, 0) = xA(1, frame) ; cell(1, 1) = xB(1, frame) ; cell(1, 2) = xC(1, frame) ;
      cell(2, 0) = xA(2, frame) ; cell(2, 1) = xB(2, frame) ; cell(2, 2) = xC(2, frame) ;
      cellInv = solve(cell) ;
    }
    for(int atom = 0; atom < natom; atom++){
      if(!Rcpp::NumericVector::is_na(xPx(atom, frame)) &&
         !Rcpp::NumericVector::is_na(yPx(atom, 0)) ){
        double dx = xPx(atom, frame) - yPx(atom, 0) ;
        double dy = xPy(atom, frame) - yPy(atom, 0) ;
        double dz = xPz(atom, frame) - yPz(atom, 0) ;
        if(usePBC[0] && (xPBC[0] || xPBC[1] || xPBC[2])){
          applyPBC2(dx, dy, dz, cell, cellInv, xPBC) ;
        }
        rmsd[frame] = dx*dx + dy*dy + dz*dz ;
        n++ ;
      }
    }
    rmsd[frame] = sqrt(rmsd[frame]/n );
  }
  
  return(rmsd) ;
}