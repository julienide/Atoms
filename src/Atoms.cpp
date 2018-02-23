#include "Atoms.h"

Atoms::Atoms()
{
}

Rcpp::S4 Atoms::asS4()
{
  setAtomicProperties() ;
  Rcpp::DataFrame bonds = Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Batm1,
    Rcpp::Named("atm2") = Batm2 ) ;
  Rcpp::DataFrame angles = Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Aatm1,
    Rcpp::Named("atm2") = Aatm2,
    Rcpp::Named("atm3") = Aatm3 ) ;
  Rcpp::DataFrame dihedrals = Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Datm1,
    Rcpp::Named("atm2") = Datm2,
    Rcpp::Named("atm3") = Datm3,
    Rcpp::Named("atm4") = Datm4 ) ;
  Rcpp::DataFrame impropers = Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Iatm1,
    Rcpp::Named("atm2") = Iatm2,
    Rcpp::Named("atm3") = Iatm3,
    Rcpp::Named("atm4") = Iatm4 ) ;
  Rcpp::NumericVector aMat = Rcpp::wrap(a) ;
  Rcpp::NumericVector bMat = Rcpp::wrap(b) ;
  Rcpp::NumericVector cMat = Rcpp::wrap(c) ;
  Rcpp::NumericVector xMat = Rcpp::wrap(x) ;
  Rcpp::NumericVector yMat = Rcpp::wrap(y) ;
  Rcpp::NumericVector zMat = Rcpp::wrap(z) ;
  Rcpp::NumericVector vxMat = Rcpp::wrap(vx) ;
  Rcpp::NumericVector vyMat = Rcpp::wrap(vy) ;
  Rcpp::NumericVector vzMat = Rcpp::wrap(vz) ;
  Rcpp::NumericVector fxMat = Rcpp::wrap(fx) ;
  Rcpp::NumericVector fyMat = Rcpp::wrap(fy) ;
  Rcpp::NumericVector fzMat = Rcpp::wrap(fz) ;
  aMat.attr("dim") = Rcpp::IntegerVector::create(3, nframe) ;
  bMat.attr("dim") = Rcpp::IntegerVector::create(3, nframe) ;
  cMat.attr("dim") = Rcpp::IntegerVector::create(3, nframe) ;
  xMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
  yMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
  zMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
  if(vx.size()){
    vxMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
    vyMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
    vzMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
  }
  if(fx.size()){
    fxMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
    fyMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
    fzMat.attr("dim") = Rcpp::IntegerVector::create(natom, nframe) ;
  }
  int current = 0 ;
  if(xMat.size())
    current = 1 ;
  Rcpp::S4 Obj("Atoms") ;
  Obj.slot("atoms") = atoms ;
  Obj.slot("bonds") = bonds ;
  Obj.slot("angles") = angles ;
  Obj.slot("dihedrals") = dihedrals ;
  Obj.slot("impropers") = impropers ;
  Obj.slot("current") = Rcpp::IntegerVector::create(current) ;
  Obj.slot("pbc") = pbc ;
  Obj.slot("a") = aMat ;
  Obj.slot("b") = bMat ;
  Obj.slot("c") = cMat ;
  Obj.slot("x") = xMat ;
  Obj.slot("y") = yMat ;
  Obj.slot("z") = zMat ;
  return Obj ;
}
