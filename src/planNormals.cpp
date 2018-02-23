#include "planNormals.h"

//' @rdname planNormals
//' @export
// [[Rcpp::export(name = ".planNormalsAtoms")]]
Rcpp::S4 planNormalsAtoms(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3)
{
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  int nframe = px.ncol() ;

  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  int nvec = atm1.size() ;
  Rcpp::NumericMatrix xout(nvec, nframe) ;
  Rcpp::NumericMatrix yout(nvec, nframe) ;
  Rcpp::NumericMatrix zout(nvec, nframe) ;
  Rcpp::NumericMatrix aout(3, nframe) ;
  Rcpp::NumericMatrix bout(3, nframe) ;
  Rcpp::NumericMatrix cout(3, nframe) ;

  for(int frame = 0; frame < nframe; frame++)
  {
    aout(0, frame) = 1.0 ;
    bout(1, frame) = 1.0 ;
    cout(2, frame) = 1.0 ;
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(int vec = 0; vec < nvec; vec++)
    {
      double ux = px(atm2[vec] - 1, frame) - px(atm1[vec] - 1, frame) ;
      double uy = py(atm2[vec] - 1, frame) - py(atm1[vec] - 1, frame) ;
      double uz = pz(atm2[vec] - 1, frame) - pz(atm1[vec] - 1, frame) ;
      applyPBC2(ux, uy, uz, cell, cellInv, pbc) ;
      double un = sqrt(ux*ux + uy*uy + uz*uz) ;
      ux = ux/un ;
      uy = uy/un ;
      uz = uz/un ;
      

      double vx = px(atm3[vec] - 1, frame) - px(atm1[vec] - 1, frame) ;
      double vy = py(atm3[vec] - 1, frame) - py(atm1[vec] - 1, frame) ;
      double vz = pz(atm3[vec] - 1, frame) - pz(atm1[vec] - 1, frame) ;
      applyPBC2(vx, vy, vz, cell, cellInv, pbc) ;
      double vn = sqrt(vx*vx + vy*vy + vz*vz) ;
      vx = vx/vn ;
      vy = vy/vn ;
      vz = vz/vn ;

      double nx = uy*vz - uz*vy ;
      double ny = uz*vx - ux*vz ;
      double nz = ux*vy - uy*vx ;
      double nn = sqrt(nx*nx + ny*ny + nz*nz) ;
      nx = nx/nn ;
      ny = ny/nn ;
      nz = nz/nn ;
      
      xout(vec, frame) = nx ;
      yout(vec, frame) = ny ;
      zout(vec, frame) = nz ;
    }
  }

  Rcpp::S4 Obj("Atoms") ;
  Rcpp::CharacterVector atmtype(nvec, "Xx") ;
  Obj.slot("atoms") = Rcpp::DataFrame::create(Rcpp::Named("atmtype") = atmtype) ;
  Obj.slot("a") = aout ;
  Obj.slot("b") = bout ;
  Obj.slot("c") = cout ;
  Obj.slot("x") = xout ;
  Obj.slot("y") = yout ;
  Obj.slot("z") = zout ;
  Obj.slot("current") = Rcpp::IntegerVector::create(1) ;
  
  return Obj ;
}