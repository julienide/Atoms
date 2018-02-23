#include "join.h"

//' @rdname join
//' @export
// [[Rcpp::export(name = ".joinAtoms")]]
void joinAtoms(
    Rcpp::S4& x,
    const Rcpp::LogicalVector& multi)
{
  bool _multi = multi[0] ;

  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  Rcpp::DataFrame bonds = x.slot("bonds") ;
  Rcpp::IntegerVector Batm1 = bonds["atm1"] ;
  Rcpp::IntegerVector Batm2 = bonds["atm2"] ;
  int nframe = px.ncol() ;
  int nbonds = Batm1.size() ;
  
  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  Progress pb(nframe, _multi) ;
  for(int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    
    for(int bond = 0; bond < nbonds; bond++){
      double dx = px[Batm2[bond] - 1] - px[Batm1[bond] - 1] ;
      double dy = py[Batm2[bond] - 1] - py[Batm1[bond] - 1] ;
      double dz = pz[Batm2[bond] - 1] - pz[Batm1[bond] - 1] ;
      applyPBC2(dx, dy, dz, cell, cellInv, pbc) ;
      px[Batm2[bond] - 1] = dx + px[Batm1[bond] - 1] ;
      py[Batm2[bond] - 1] = dy + py[Batm1[bond] - 1] ;
      pz[Batm2[bond] - 1] = dz + pz[Batm1[bond] - 1] ;
    }
    
    pb.increment() ;
  }
}