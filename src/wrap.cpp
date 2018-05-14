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
      double cxf = cx[res]*cellInv(0, 0) + cy[res]*cellInv(0, 1) + cz[res]*cellInv(0, 2) ;
      double cyf = cx[res]*cellInv(1, 0) + cy[res]*cellInv(1, 1) + cz[res]*cellInv(1, 2) ;
      double czf = cx[res]*cellInv(2, 0) + cy[res]*cellInv(2, 1) + cz[res]*cellInv(2, 2) ;
      if(pbc[0])
      {
        cxf = -std::floor(cxf) ;
      } else {
        cxf = 0.0 ;
      }
      if(pbc[1])
      {
        cyf = -std::floor(cyf) ;
      } else {
        cyf = 0.0 ;
      }
      if(pbc[2])
      {
        czf = -std::floor(czf) ;
      } else {
        czf = 0.0 ;
      }
      cx[res] = cxf*cell(0, 0) + cyf*cell(0, 1) + czf*cell(0, 2) ;
      cy[res] = cxf*cell(1, 0) + cyf*cell(1, 1) + czf*cell(1, 2) ;
      cz[res] = cxf*cell(2, 0) + cyf*cell(2, 1) + czf*cell(2, 2) ;
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
