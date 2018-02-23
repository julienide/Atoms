#include "freeVolume.h"

const double pi = std::atan(1.0)*4 ;

double VSphere(double r){
  double V = r*r*r*pi*4.0/3.0 ;
  return V ;
}

double HLens(double R1, double R2, double d)
{
  double h = (R2 - R1 + d)*(R2 + R1 - d)/(2.0*d) ;
  return h ;
}

double VLens(double r, double h)
{
  double V = h*h*(3.0*r - h)*pi/3.0 ;
  return V ;
}

//' @rdname freeVolume
//' @export
// [[Rcpp::export(name = "freeVolume.Atoms")]]
Rcpp::NumericVector freeVolumeAtoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable< Rcpp::NumericVector > radius = R_NilValue)
{
  std::vector< double > r ;
  if(radius.isNull()){
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    if(atmprop.containsElementNamed("radius")){
      r = Rcpp::as< std::vector< double > >(atmprop["radius"]) ;
    } else {
      Rcpp::Environment PeriodicTableEnv("package:PeriodicTable") ;
      Rcpp::Function rvdw = PeriodicTableEnv["rvdw"] ;
      r = Rcpp::as< std::vector< double > >(rvdw(x)) ;
    }
  } else {
    r = Rcpp::as< std::vector< double > >(radius) ;
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    int rsize = r.size() ;
    if(rsize != atmprop.nrows()){
      Rcpp::stop("'radius' must be of length ", atmprop.nrows()) ;
    }
  }
  
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  int nframe = px.ncol() ;
  
  Rcpp::Environment BaseEnv("package:base") ;
  Rcpp::Function solve = BaseEnv["solve"] ;

  Rcpp::NumericVector freev(nframe) ;
  
  Progress pb(nframe, true) ;
  for(int frame = 0; frame != nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(int at1 = 0; at1 != px.nrow(); at1++) {
      freev[frame] = freev[frame] + VSphere(r[at1]) ;
      for(int at2 = at1 + 1; at2 != px.nrow(); at2++) {
        double dx = px(at2, frame) - px(at1, frame) ;
        double dy = py(at2, frame) - py(at1, frame) ;
        double dz = pz(at2, frame) - pz(at1, frame) ;
        if(pbc[0] || pbc[1] || pbc[2])
        {
          applyPBC2(dx, dy, dz, cell, cellInv, pbc) ;
        }
        double d = sqrt(dx*dx + dy*dy + dz*dz) ;
        double dbond = r[at1] + r[at2] ;
        if(d < dbond)
        {
          if(d > abs(r[at1] - r[at2]))
          {
            // Intersecting spheres
            double h1 = HLens(r[at1], r[at2], d) ;
            double h2 = HLens(r[at2], r[at1], d) ;
            freev[frame] = freev[frame] - VLens(r[at1], h1) - VLens(r[at2], h2) ;
          } else {
            // at1 completly inside at2. This should never append.
            // if(r[at1] < r[at2])
            // {
            //   // The volume of at1 must be substracted only for the first completly overlapping neighbor!
            //   freev[frame] = freev[frame] - VSphere(r[at1]) ;
            // }
          }
          
          // double h1 ;
          // if(d >= r[at1]){
          //   h1 = (dbond - d)/2.0 ;
          // } else {
          //   h1 = r[at1] - d/2.0 ;
          // }
          // double h2 ;
          // if(d >= r[at2]){
          //   h2 = (dbond - d)/2.0 ;
          // } else {
          //   h2 = r[at2] - d/2.0 ;
          // }
          // freev[frame] = freev[frame] - VLens(r[at1], h1) - VLens(r[at2], h2) ;
        }
      }
    }
    double V =
      cell(0, 0)*(cell(1, 1)*cell(2, 2) - cell(2, 1)*cell(1, 2)) +
      cell(0, 1)*(cell(2, 1)*cell(0, 2) - cell(0, 1)*cell(2, 2)) +
      cell(0, 2)*(cell(0, 1)*cell(1, 2) - cell(1, 1)*cell(0, 2)) ;
    if(freev[frame] < 0){
      freev[frame] = 0.0 ;
    } else {
      freev[frame] = V - freev[frame] ;
    }
    pb.increment() ;
  }
  
  return freev ;
}

//' @rdname freeVolume
//' @export
// [[Rcpp::export(name = "freeVolume2.Atoms")]]
Rcpp::NumericVector freeVolume2Atoms(
    const Rcpp::S4& x,
    const Rcpp::Nullable< Rcpp::NumericVector > radius = R_NilValue)
{
  std::vector< double > r ;
  if(radius.isNull()){
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    if(atmprop.containsElementNamed("radius")){
      r = Rcpp::as< std::vector< double > >(atmprop["radius"]) ;
    } else {
      Rcpp::Environment PeriodicTableEnv("package:PeriodicTable") ;
      Rcpp::Function rvdw = PeriodicTableEnv["rvdw"] ;
      r = Rcpp::as< std::vector< double > >(rvdw(x)) ;
    }
  } else {
    r = Rcpp::as< std::vector< double > >(radius) ;
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    int rsize = r.size() ;
    if(rsize != atmprop.nrows()){
      Rcpp::stop("'radius' must be of length ", atmprop.nrows()) ;
    }
  }

  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  int nframe = px.ncol() ;

  Rcpp::Environment BaseEnv("package:base") ;
  Rcpp::Function solve = BaseEnv["solve"] ;

  Rcpp::NumericVector freev(nframe) ;

  Progress pb(nframe, true) ;
  for(int frame = 0; frame != nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(int at1 = 0; at1 != px.nrow(); at1++) {
      freev[frame] = freev[frame] + VSphere(r[at1]) ;
      for(int at2 = at1 + 1; at2 != px.nrow(); at2++) {
        double dx = px(at2, frame) - px(at1, frame) ;
        double dy = py(at2, frame) - py(at1, frame) ;
        double dz = pz(at2, frame) - pz(at1, frame) ;
        if(pbc[0] || pbc[1] || pbc[2])
        {
          applyPBC2(dx, dy, dz, cell, cellInv, pbc) ;
        }
        double d = sqrt(dx*dx + dy*dy + dz*dz) ;
        double dbond = r[at1] + r[at2] ;
        if(d < dbond)
        {
          double h ;
          if(d >= r[at1]){
            h = (dbond - d)/2.0 ;
          } else {
            h = r[at1] - d/2.0 ;
          }
          freev[frame] = freev[frame] - VLens(r[at1], h) ;
          if(d >= r[at2]){
            h = (dbond - d)/2.0 ;
          } else {
            h = r[at2] - d/2.0 ;
          }
          freev[frame] = freev[frame] - VLens(r[at2], h) ;
        }
      }
    }
    double V =
      cell(0, 0)*(cell(1, 1)*cell(2, 2) - cell(2, 1)*cell(1, 2)) +
      cell(0, 1)*(cell(2, 1)*cell(0, 2) - cell(0, 1)*cell(2, 2)) +
      cell(0, 2)*(cell(0, 1)*cell(1, 2) - cell(1, 1)*cell(0, 2)) ;
    if(freev[frame] < 0.0){
      freev[frame] = 0.0 ;
    } else {
      freev[frame] = V - freev[frame] ;
    }
    pb.increment() ;
  }

  return freev ;
}
