#include <Rcpp.h>
#include "com.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

//' @rdname com
//' @export
// [[Rcpp::export(name = ".comAtoms")]]
Rcpp::S4 comAtoms(
    const Rcpp::S4& x,
    const Rcpp::NumericVector& mass,
    const Rcpp::IntegerVector& factor){
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  Rcpp::CharacterVector lvls = factor.attr("levels") ;

  int natom = px.nrow() ;
  int nframe = px.ncol() ;
  int nres = lvls.size() ;

  Rcpp::NumericVector mcom(nres, 0.0) ;
  Rcpp::NumericMatrix xcom(nres, nframe) ;
  Rcpp::NumericMatrix ycom(nres, nframe) ;
  Rcpp::NumericMatrix zcom(nres, nframe) ;

  for(int atom = 0; atom < natom; atom++)
  {
    int res = factor(atom) ;
    if(res != NA_INTEGER){
      res = res - 1 ;
      mcom(res) = mcom(res) + mass(atom) ;
    }
  }

  Progress pb(natom, true) ;
  for(int atom = 0; atom < natom; atom++)
  {
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    int res = factor(atom) ;
    if(res != NA_INTEGER){
      res = res - 1 ;
      for(int frame = 0; frame < nframe; frame++)
      {
        xcom(res, frame) = xcom(res, frame) + px(atom, frame)*mass(atom)/mcom(res) ;
        ycom(res, frame) = ycom(res, frame) + py(atom, frame)*mass(atom)/mcom(res) ;
        zcom(res, frame) = zcom(res, frame) + pz(atom, frame)*mass(atom)/mcom(res) ;
      }
    }
    pb.increment() ;
  }

  Rcpp::S4 Obj ("Atoms") ;
  Obj.slot("a") = x.slot("a") ;
  Obj.slot("b") = x.slot("b") ;
  Obj.slot("c") = x.slot("c") ;
  Obj.slot("pbc") = x.slot("pbc") ;
  Obj.slot("current") = Rcpp::IntegerVector::create(1) ;
  Obj.slot("x") = xcom ;
  Obj.slot("y") = ycom ;
  Obj.slot("z") = zcom ;
  Obj.slot("atoms") = Rcpp::DataFrame::create(
    Rcpp::Named("atmname") = "Xx",
    Rcpp::Named("atmtype") = "Xx",
    Rcpp::Named("resname") = lvls,
    Rcpp::Named("mass") = mcom,
    Rcpp::Named("stringsAsFactors") = false) ;

  return Obj ;
}
