#include "vectOrientation.h"

const double pi = std::atan(1.0)*4 ;

double norm(
    const double& x, const double& y, const double& z){
  double n = sqrt(x*x + y*y + z*z) ;
  return n ;
}

void normalize(
    double& x, double& y, double& z){
  double n = norm(x, y, z) ;
  x = x/n ;
  y = y/n ;
  z = z/n ;
}

std::vector< double > distVect(
    const double& x1, const double& y1, const double& z1,
    const double& x2, const double& y2, const double& z2,
    const Rcpp::NumericMatrix& cell, const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc){
  std::vector< double > d(3) ;
  d[0] = x2 - x1 ;
  d[1] = y2 - y1 ;
  d[2] = z2 - z1 ;
  applyPBC2(d[0], d[1], d[2], cell, cellInv, pbc) ;
  normalize(d[0], d[1], d[2]) ;
  return d ;
}

std::vector< double > planNormal(
    const double& x1, const double& y1, const double& z1,
    const double& x2, const double& y2, const double& z2,
    const double& x3, const double& y3, const double& z3,
    const Rcpp::NumericMatrix& cell, const Rcpp::NumericMatrix& cellInv,
    const Rcpp::LogicalVector& pbc){
  std::vector< double > u = distVect(x1, y1, z1, x2, y2, z2, cell, cellInv, pbc) ;
  std::vector< double > v = distVect(x1, y1, z1, x3, y3, z3, cell, cellInv, pbc) ;
  std::vector< double > n(3) ;
  n[0] = u[1]*v[2] - u[2]*v[1] ;
  n[1] = u[2]*v[0] - u[0]*v[2] ;
  n[2] = u[0]*v[1] - u[1]*v[0] ;
  normalize(n[0], n[1], n[2]) ;
  return n ;
}

double angle(double V1X, double V1Y, double V1Z,
             double V2X, double V2Y, double V2Z){
  double A = V1X*V2X + V1Y*V2Y + V1Z*V2Z ;
  A = acos(A)*180.0/pi ;
  return A ;
}

//' @rdname vectOrientation
//' @noRd
//' @export
// [[Rcpp::export(name = ".vectOrientation")]]
Rcpp::NumericVector vectOrientation(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1){
    Rcpp::NumericMatrix px = x.slot("x") ;
    Rcpp::NumericMatrix py = x.slot("y") ;
    Rcpp::NumericMatrix pz = x.slot("z") ;
    Rcpp::NumericMatrix a = x.slot("a") ;
    Rcpp::NumericMatrix b = x.slot("b") ;
    Rcpp::NumericMatrix c = x.slot("c") ;
    Rcpp::LogicalVector pbc = x.slot("pbc") ;
    int nframe = px.ncol() ;

    Rcpp::IntegerVector V1Atm1 = V1["atm1"] ;
    Rcpp::IntegerVector V1Atm2 = V1["atm2"] ;
    // Rcpp::IntegerVector V2Atm1 = V2["atm1"] ;
    // Rcpp::IntegerVector V2Atm2 = V2["atm2"] ;
    int NV1 = V1Atm1.size() ;
    // int NV2 = V2Atm1.size() ;
    
    bool V1VectType = true ;
    Rcpp::IntegerVector V1Atm3 ;
    if(V1.containsElementNamed("atm3")){
      V1Atm3 = V1["atm3"] ;
      V1VectType = false ;
    }
// 
//     bool V2VectType = true ;
//     Rcpp::IntegerVector V2Atm3 ;
//     if(V2.containsElementNamed("atm3")){
//       V2Atm3 = V2["atm3"] ;
//       V2VectType = false ;
//     }
    
    Rcpp::Environment baseEnv("package:base") ;
    Rcpp::Function solve = baseEnv["solve"] ;
    
    std::vector< double > AAll ;
    
    for(int frame = 0; frame < nframe; frame++)
    {
      Rcpp::NumericMatrix cell(3, 3) ;
      cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
      cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
      cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
      Rcpp::NumericMatrix cellInv = solve(cell) ;
      for(int V1Row = 0; V1Row < NV1; V1Row++){
        std::vector< double > V1Vect(3) ;
        if(V1VectType){
          V1Vect = distVect(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                            px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1], 
                            cell, cellInv, pbc) ;
        } else {
          V1Vect = planNormal(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                              px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1],
                              px[V1Atm3[V1Row] - 1], py[V1Atm3[V1Row] - 1], pz[V1Atm3[V1Row] - 1],
                              cell, cellInv, pbc) ;
        }
        for(int V2Row = V1Row + 1; V2Row < NV1; V2Row++){
          std::vector< double > V2Vect(3) ;
          if(V1VectType){
            V2Vect = distVect(px[V1Atm1[V2Row] - 1], py[V1Atm1[V2Row] - 1], pz[V1Atm1[V2Row] - 1],
                              px[V1Atm2[V2Row] - 1], py[V1Atm2[V2Row] - 1], pz[V1Atm2[V2Row] - 1], 
                              cell, cellInv, pbc) ;
          } else {
            V2Vect = planNormal(px[V1Atm1[V2Row] - 1], py[V1Atm1[V2Row] - 1], pz[V1Atm1[V2Row] - 1],
                                px[V1Atm2[V2Row] - 1], py[V1Atm2[V2Row] - 1], pz[V1Atm2[V2Row] - 1],
                                px[V1Atm3[V2Row] - 1], py[V1Atm3[V2Row] - 1], pz[V1Atm3[V2Row] - 1],
                                cell, cellInv, pbc) ;
          }
          double A = angle(V1Vect[0], V1Vect[1], V1Vect[2],
                           V2Vect[0], V2Vect[1], V2Vect[2]) ;
          AAll.push_back(A) ;
        }
      }
    }
    Rcpp::NumericVector AAllOut = Rcpp::wrap(AAll) ;
    AAllOut.attr("dim") = Rcpp::IntegerVector::create(NV1, NV1, nframe) ;
    return AAllOut ;
}

//' @rdname vectOrientation
//' @noRd
//' @export
// [[Rcpp::export(name = ".vectOrientationAxis")]]
Rcpp::NumericVector vectOrientationAxis(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::NumericVector& V2){
    Rcpp::NumericMatrix px = x.slot("x") ;
    Rcpp::NumericMatrix py = x.slot("y") ;
    Rcpp::NumericMatrix pz = x.slot("z") ;
    Rcpp::NumericMatrix a = x.slot("a") ;
    Rcpp::NumericMatrix b = x.slot("b") ;
    Rcpp::NumericMatrix c = x.slot("c") ;
    Rcpp::LogicalVector pbc = x.slot("pbc") ;
    int nframe = px.ncol() ;
    
    Rcpp::IntegerVector V1Atm1 = V1["atm1"] ;
    Rcpp::IntegerVector V1Atm2 = V1["atm2"] ;
    int NV1 = V1Atm1.size() ;

    bool V1VectType = true ;
    Rcpp::IntegerVector V1Atm3 ;
    if(V1.containsElementNamed("atm3")){
      V1Atm3 = V1["atm3"] ;
      V1VectType = false ;
    }

    std::vector< double > V2Vect(3) ;
    V2Vect[0] = V2[0] ;
    V2Vect[1] = V2[1] ;
    V2Vect[2] = V2[2] ;
    normalize(V2Vect[0], V2Vect[1], V2Vect[2]) ;
    
    Rcpp::Environment baseEnv("package:base") ;
    Rcpp::Function solve = baseEnv["solve"] ;
    
    std::vector< double > AAll ;
    
    for(int frame = 0; frame < nframe; frame++)
    {
      Rcpp::NumericMatrix cell(3, 3) ;
      cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
      cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
      cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
      Rcpp::NumericMatrix cellInv = solve(cell) ;
      for(int V1Row = 0; V1Row < NV1; V1Row++){
        std::vector< double > V1Vect(3) ;
        if(V1VectType){
          V1Vect = distVect(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                            px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1], 
                                                                            cell, cellInv, pbc) ;
        } else {
          V1Vect = planNormal(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                              px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1],
                                                                              px[V1Atm3[V1Row] - 1], py[V1Atm3[V1Row] - 1], pz[V1Atm3[V1Row] - 1],
                                                                                                                              cell, cellInv, pbc) ;
        }
        double A = angle(V1Vect[0], V1Vect[1], V1Vect[2],
                         V2Vect[0], V2Vect[1], V2Vect[2]) ;
        AAll.push_back(A) ;
      }
    }
    Rcpp::NumericVector AAllOut = Rcpp::wrap(AAll) ;
    AAllOut.attr("dim") = Rcpp::IntegerVector::create(NV1, nframe) ;
    return AAllOut ;
}

//' @rdname vectOrientation
//' @noRd
//' @export
// [[Rcpp::export(name = ".vectOrientationPairwise")]]
Rcpp::NumericVector vectOrientationPairwise(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::DataFrame& V2){
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  int nframe = px.ncol() ;
  
  Rcpp::IntegerVector V1Atm1 = V1["atm1"] ;
  Rcpp::IntegerVector V1Atm2 = V1["atm2"] ;
  Rcpp::IntegerVector V2Atm1 = V2["atm1"] ;
  Rcpp::IntegerVector V2Atm2 = V2["atm2"] ;
  int NV1 = V1Atm1.size() ;

  bool V1VectType = true ;
  Rcpp::IntegerVector V1Atm3 ;
  if(V1.containsElementNamed("atm3")){
    V1Atm3 = V1["atm3"] ;
    V1VectType = false ;
  }
  
  bool V2VectType = true ;
  Rcpp::IntegerVector V2Atm3 ;
  if(V2.containsElementNamed("atm3")){
    V2Atm3 = V2["atm3"] ;
    V2VectType = false ;
  }
  
  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  std::vector< double > AAll ;
  
  for(int frame = 0; frame < nframe; frame++)
  {
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(int VRow = 0; VRow < NV1; VRow++){
      std::vector< double > V1Vect(3) ;
      if(V1VectType){
        V1Vect = distVect(px[V1Atm1[VRow] - 1], py[V1Atm1[VRow] - 1], pz[V1Atm1[VRow] - 1],
                          px[V1Atm2[VRow] - 1], py[V1Atm2[VRow] - 1], pz[V1Atm2[VRow] - 1], 
                          cell, cellInv, pbc) ;
      } else {
        V1Vect = planNormal(px[V1Atm1[VRow] - 1], py[V1Atm1[VRow] - 1], pz[V1Atm1[VRow] - 1],
                            px[V1Atm2[VRow] - 1], py[V1Atm2[VRow] - 1], pz[V1Atm2[VRow] - 1],
                            px[V1Atm3[VRow] - 1], py[V1Atm3[VRow] - 1], pz[V1Atm3[VRow] - 1],
                            cell, cellInv, pbc) ;
      }
      std::vector< double > V2Vect(3) ;
      if(V2VectType){
        V2Vect = distVect(px[V2Atm1[VRow] - 1], py[V2Atm1[VRow] - 1], pz[V2Atm1[VRow] - 1],
                          px[V2Atm2[VRow] - 1], py[V2Atm2[VRow] - 1], pz[V2Atm2[VRow] - 1],
                          cell, cellInv, pbc) ;
      } else {
        V2Vect = planNormal(px[V2Atm1[VRow] - 1], py[V2Atm1[VRow] - 1], pz[V2Atm1[VRow] - 1],
                            px[V2Atm2[VRow] - 1], py[V2Atm2[VRow] - 1], pz[V2Atm2[VRow] - 1],
                            px[V2Atm3[VRow] - 1], py[V2Atm3[VRow] - 1], pz[V2Atm3[VRow] - 1],
                            cell, cellInv, pbc) ;
      }
      double A = angle(V1Vect[0], V1Vect[1], V1Vect[2],
                       V2Vect[0], V2Vect[1], V2Vect[2]) ;
      AAll.push_back(A) ;
    }
  }
  Rcpp::NumericVector AAllOut = Rcpp::wrap(AAll) ;
  AAllOut.attr("dim") = Rcpp::IntegerVector::create(NV1, nframe) ;
  return AAllOut ;
}

//' @rdname vectOrientation
//' @noRd
//' @export
// [[Rcpp::export(name = ".vectOrientationNotPairwise")]]
Rcpp::NumericVector vectOrientationNotPairwise(
    const Rcpp::S4& x,
    const Rcpp::DataFrame& V1,
    const Rcpp::DataFrame& V2){
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  int nframe = px.ncol() ;
  
  Rcpp::IntegerVector V1Atm1 = V1["atm1"] ;
  Rcpp::IntegerVector V1Atm2 = V1["atm2"] ;
  Rcpp::IntegerVector V2Atm1 = V2["atm1"] ;
  Rcpp::IntegerVector V2Atm2 = V2["atm2"] ;
  int NV1 = V1Atm1.size() ;
  int NV2 = V2Atm1.size() ;
  
  bool V1VectType = true ;
  Rcpp::IntegerVector V1Atm3 ;
  if(V1.containsElementNamed("atm3")){
    V1Atm3 = V1["atm3"] ;
    V1VectType = false ;
  }
  
  bool V2VectType = true ;
  Rcpp::IntegerVector V2Atm3 ;
  if(V2.containsElementNamed("atm3")){
    V2Atm3 = V2["atm3"] ;
    V2VectType = false ;
  }
  
  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;
  
  std::vector< double > AAll ;
  
  for(int frame = 0; frame < nframe; frame++)
  {
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    for(int V1Row = 0; V1Row < NV1; V1Row++){
      std::vector< double > V1Vect(3) ;
      if(V1VectType){
        V1Vect = distVect(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                          px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1], 
                                                                          cell, cellInv, pbc) ;
      } else {
        V1Vect = planNormal(px[V1Atm1[V1Row] - 1], py[V1Atm1[V1Row] - 1], pz[V1Atm1[V1Row] - 1],
                            px[V1Atm2[V1Row] - 1], py[V1Atm2[V1Row] - 1], pz[V1Atm2[V1Row] - 1],
                                                                            px[V1Atm3[V1Row] - 1], py[V1Atm3[V1Row] - 1], pz[V1Atm3[V1Row] - 1],
                                                                                                                            cell, cellInv, pbc) ;
      }
      for(int V2Row = 0; V2Row < NV1; V2Row++){
        std::vector< double > V2Vect(3) ;
        if(V2VectType){
          V2Vect = distVect(px[V2Atm1[V2Row] - 1], py[V2Atm1[V2Row] - 1], pz[V2Atm1[V2Row] - 1],
                            px[V2Atm2[V2Row] - 1], py[V2Atm2[V2Row] - 1], pz[V2Atm2[V2Row] - 1], 
                            cell, cellInv, pbc) ;
        } else {
          V2Vect = planNormal(px[V2Atm1[V2Row] - 1], py[V2Atm1[V2Row] - 1], pz[V2Atm1[V2Row] - 1],
                              px[V2Atm2[V2Row] - 1], py[V2Atm2[V2Row] - 1], pz[V2Atm2[V2Row] - 1],
                              px[V2Atm3[V2Row] - 1], py[V2Atm3[V2Row] - 1], pz[V2Atm3[V2Row] - 1],
                              cell, cellInv, pbc) ;
        }
        double A = angle(V1Vect[0], V1Vect[1], V1Vect[2],
                         V2Vect[0], V2Vect[1], V2Vect[2]) ;
        AAll.push_back(A) ;
      }
    }
  }
  Rcpp::NumericVector AAllOut = Rcpp::wrap(AAll) ;
  AAllOut.attr("dim") = Rcpp::IntegerVector::create(NV1, NV2, nframe) ;
  return AAllOut ;
}
