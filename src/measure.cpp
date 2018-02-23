#include "measure.h"

//' @rdname measure
//' @noRd
//' @export
// [[Rcpp::export(name = ".measureBonds")]]
Rcpp::NumericMatrix measureBonds(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2){
  if(atm1.size() != atm2.size()){
    Rcpp::stop("'atm1' and 'atm2' must have the same length") ;
  }
  
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;

  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  unsigned int nbond = atm1.size() ;
  unsigned int nframe = px.ncol() ;
  
  Progress pb(nframe, true) ;
  
  Rcpp::NumericMatrix B(nbond, nframe) ;
  for(unsigned int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(unsigned int bond = 0; bond < nbond; bond++){
      double dx = px(atm2[bond] - 1, frame) - px(atm1[bond] - 1, frame) ;
      double dy = py(atm2[bond] - 1, frame) - py(atm1[bond] - 1, frame) ;
      double dz = pz(atm2[bond] - 1, frame) - pz(atm1[bond] - 1, frame) ;
      applyPBC2(dx, dy, dz, cell, cellInv, pbc) ;
      B(bond, frame) = sqrt(dx*dx + dy*dy + dz*dz) ;
    }
    pb.increment() ;
  }
  return B ;
}

//' @rdname measure
//' @noRd
//' @export
// [[Rcpp::export(name = ".measureAngles")]]
Rcpp::NumericMatrix measureAngles(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3){
  if(atm1.size() != atm2.size() ||
     atm1.size() != atm3.size() ){
    Rcpp::stop("'atm1', 'atm2' and 'atm3' must have the same length") ;
  }
  
  const double pi = std::atan(1.0)*4 ;
  
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  
  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  unsigned int nangle = atm1.size() ;
  unsigned int nframe = px.ncol() ;
  
  Progress pb(nframe, true) ;
  
  Rcpp::NumericMatrix A(nangle, nframe) ;
  for(unsigned int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(unsigned int angle = 0; angle < nangle; angle++){
      double d21x = px(atm1[angle] - 1, frame) - px(atm2[angle] - 1, frame) ;
      double d21y = py(atm1[angle] - 1, frame) - py(atm2[angle] - 1, frame) ;
      double d21z = pz(atm1[angle] - 1, frame) - pz(atm2[angle] - 1, frame) ;
      double d32x = px(atm3[angle] - 1, frame) - px(atm2[angle] - 1, frame) ;
      double d32y = py(atm3[angle] - 1, frame) - py(atm2[angle] - 1, frame) ;
      double d32z = pz(atm3[angle] - 1, frame) - pz(atm2[angle] - 1, frame) ;
      applyPBC2(d21x, d21y, d21z, cell, cellInv, pbc) ;
      applyPBC2(d32x, d32y, d32z, cell, cellInv, pbc) ;
      A(angle, frame) = d21x*d32x + d21y*d32y + d21z*d32z ;
      A(angle, frame) = A(angle, frame)/
        sqrt(d21x*d21x + d21y*d21y + d21z*d21z) ;
      A(angle, frame) = A(angle, frame)/
        sqrt(d32x*d32x + d32y*d32y + d32z*d32z) ;
      A(angle, frame) = std::acos(A(angle, frame))*180.0/pi ;
    }
    pb.increment() ;
  }
  return A ;
}

//' @rdname measure
//' @noRd
//' @export
// [[Rcpp::export(name = ".measureDihedrals")]]
Rcpp::NumericMatrix measureDihedrals(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& atm1,
    const Rcpp::IntegerVector& atm2,
    const Rcpp::IntegerVector& atm3,
    const Rcpp::IntegerVector& atm4){
  if(atm1.size() != atm2.size() ||
     atm1.size() != atm3.size() ||
     atm1.size() != atm4.size()){
    Rcpp::stop("'atm1', 'atm2', 'atm3' and 'atm4' must have the same length") ;
  }
  
  const double pi = std::atan(1.0)*4 ;

  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;

  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  
  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  unsigned int ndihedrals = atm1.size() ;
  unsigned int nframe = px.ncol() ;

  Progress pb(nframe, true) ;
  
  Rcpp::NumericMatrix D(ndihedrals, nframe) ;
  for(unsigned int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(unsigned int dihedral = 0; dihedral < ndihedrals; dihedral++){
      double ux = px[atm2[dihedral] - 1] - px[atm1[dihedral] - 1] ;
      double uy = py[atm2[dihedral] - 1] - py[atm1[dihedral] - 1] ;
      double uz = pz[atm2[dihedral] - 1] - pz[atm1[dihedral] - 1] ;
      
      double vx = px[atm3[dihedral] - 1] - px[atm2[dihedral] - 1] ;
      double vy = py[atm3[dihedral] - 1] - py[atm2[dihedral] - 1] ;
      double vz = pz[atm3[dihedral] - 1] - pz[atm2[dihedral] - 1] ;
      
      double wx = px[atm4[dihedral] - 1] - px[atm3[dihedral] - 1] ;
      double wy = py[atm4[dihedral] - 1] - py[atm3[dihedral] - 1] ;
      double wz = pz[atm4[dihedral] - 1] - pz[atm3[dihedral] - 1] ;
      
      applyPBC2(ux, uy, uz, cell, cellInv, pbc) ;
      applyPBC2(vx, vy, vz, cell, cellInv, pbc) ;
      applyPBC2(wx, wy, wz, cell, cellInv, pbc) ;
      
      double un = sqrt(ux*ux + uy*uy + uz*uz) ;
      double vn = sqrt(vx*vx + vy*vy + vz*vz) ;
      double wn = sqrt(wx*wx + wy*wy + wz*wz) ;
      
      ux = ux/un ;
      uy = uy/un ;
      uz = uz/un ;
      
      vx = vx/vn ;
      vy = vy/vn ;
      vz = vz/vn ;
      
      wx = wx/wn ;
      wy = wy/wn ;
      wz = wz/wn ;

      double n1x = uz*vy - uy*vz ;
      double n1y = ux*vz - uz*vx ;
      double n1z = uy*vx - ux*vy ;
      
      double n2x = vz*wy - vy*wz ;
      double n2y = vx*wz - vz*wx ;
      double n2z = vy*wx - vx*wy ;
      
      double m1x = n1z*vy - n1y*vz ;
      double m1y = n1x*vz - n1z*vx ;
      double m1z = n1y*vx - n1x*vy ;
      
      double n1n = sqrt(n1x*n1x + n1y*n1y + n1z*n1z) ;
      double n2n = sqrt(n2x*n2x + n2y*n2y + n2z*n2z) ;
      double m1n = sqrt(m1x*m1x + m1y*m1y + m1z*m1z) ;
      
      n1x = n1x/n1n ;
      n1y = n1y/n1n ;
      n1z = n1z/n1n ;
      
      n2x = n2x/n2n ;
      n2y = n2y/n2n ;
      n2z = n2z/n2n ;
      
      m1x = m1x/m1n ;
      m1y = m1y/m1n ;
      m1z = m1z/m1n ;
      
      double x = n1x*n2x + n1y*n2y + n1z*n2z ;
      double y = m1x*n2x + m1y*n2y + m1z*n2z ;
      
      D(dihedral, frame) = std::atan2(y, x)*180.0/pi ;
    }
    pb.increment() ;
  }
  return D ;
}
