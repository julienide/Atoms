#include "wrap.h"

//' @rdname wrap
//' @export
// [[Rcpp::export(name = ".wrapAtoms")]]
void wrapAtoms(
    Rcpp::S4& x,
    const Rcpp::IntegerVector& resnumb,
    const Rcpp::LogicalVector& multi)
{
  bool _multi = multi[0] ;
  
  Rcpp::CharacterVector lvls = resnumb.attr("levels") ;
  int nlvls = lvls.size() ;
  
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  int nframe = px.ncol() ;
  int natom = px.nrow() ;
  
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
    
    std::vector< double > cx(nlvls, 0.0) ;
    std::vector< double > cy(nlvls, 0.0) ;
    std::vector< double > cz(nlvls, 0.0) ;
    std::vector< int > cn(nlvls, 0) ;
    for(int atom = 0; atom < natom; atom++)
    {
      int res = resnumb[atom] - 1L ;
      cx[res] = cx[res] + px(atom, frame) ;
      cy[res] = cy[res] + py(atom, frame) ;
      cz[res] = cz[res] + pz(atom, frame) ;
      cn[res] = cn[res] + 1 ;
    }
    for(int res = 0; res < nlvls; res++)
    {
      cx[res] = cx[res]/cn[res] ;
      cy[res] = cy[res]/cn[res] ;
      cz[res] = cz[res]/cn[res] ;
      cx[res] = cx[res]*cellInv(0, 0) + cy[res]*cellInv(0, 1) + cz[res]*cellInv(0, 2) ;
      cy[res] = cx[res]*cellInv(1, 0) + cy[res]*cellInv(1, 1) + cz[res]*cellInv(1, 2) ;
      cz[res] = cx[res]*cellInv(2, 0) + cy[res]*cellInv(2, 1) + cz[res]*cellInv(2, 2) ;
      if(pbc[0])
      {
        cx[res] = -std::floor(cx[res]) ;
      } else {
        cx[res] = 0.0 ;
      }
      if(pbc[1])
      {
        cy[res] = -std::floor(cy[res]) ;
      } else {
        cy[res] = 0.0 ;
      }
      if(pbc[2])
      {
        cz[res] = -std::floor(cz[res]) ;
      } else {
        cz[res] = 0.0 ;
      }
      cx[res] = cx[res]*cell(0, 0) + cy[res]*cell(0, 1) + cz[res]*cell(0, 2) ;
      cy[res] = cx[res]*cell(1, 0) + cy[res]*cell(1, 1) + cz[res]*cell(1, 2) ;
      cz[res] = cx[res]*cell(2, 0) + cy[res]*cell(2, 1) + cz[res]*cell(2, 2) ;
    }
    for(int atom = 0; atom < natom; atom++)
    {
      int res = resnumb[atom] - 1L ;
      px(atom, frame) = px(atom, frame) + cx[res] ;
      py(atom, frame) = py(atom, frame) + cy[res] ;
      pz(atom, frame) = pz(atom, frame) + cz[res] ;
    }
    pb.increment() ;
  }
}
