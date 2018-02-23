#include <Rcpp.h>
#include "guessBonds.h"

//' @param radius a numeric vector. Atomic radii used to defined bonded atoms.
//' @param safety a numeric value used to extend the atomic radii.
//' @param usePBC a logical value. Whether to use periodic boundary conditions to compute inter-atomic distances.
//' @rdname guessTopology
//' @export
// [[Rcpp::export(name = "guessBonds.Atoms")]]
Rcpp::DataFrame guessBondsAtoms( // Check out VMD: BondSearch.C
  const Rcpp::S4& x,
  const Rcpp::Nullable< Rcpp::NumericVector > radius = R_NilValue,
  const Rcpp::NumericVector& safety = 1.2,
  const Rcpp::LogicalVector& usePBC = true){
  std::vector< double > r ;
  if(radius.isNull()){
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    if(atmprop.containsElementNamed("radius")){
      r = Rcpp::as< std::vector< double > >(atmprop["radius"]) ;
    } else {
      Rcpp::Environment PeriodicTableEnv("package:PeriodicTable") ;
      Rcpp::Function rcov = PeriodicTableEnv["rcov"] ;
      r = Rcpp::as< std::vector< double > >(rcov(x)) ;
    }
  } else {
    r = Rcpp::as< std::vector< double > >(radius) ;
    Rcpp::DataFrame atmprop = x.slot("atoms") ;
    int rsize = r.size() ;
    if(rsize != atmprop.nrows()){
      Rcpp::stop("'radius' must be of length ", atmprop.nrows()) ;
    }
  }
  if(safety.size() != 1)
    Rcpp::stop("'safety' must be a numeric vector of length 1") ;
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::IntegerVector current = x.slot("current") ;
  int frame = current[0] - 1 ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;

  Rcpp::Environment AtomsEnv("package:Atoms") ;
  Rcpp::Function getCell = AtomsEnv["cell"] ;
  Rcpp::NumericMatrix cell = getCell(x) ;

  Rcpp::Environment BaseEnv("package:base") ;
  Rcpp::Function solve = BaseEnv["solve"] ;
  Rcpp::NumericMatrix cellInv = solve(cell) ;

  std::vector<int> atm1, atm2, h, k, l ;
  for(int at1 = 0; at1 != px.nrow(); at1++) {
    for(int at2 = at1 + 1; at2 != px.nrow(); at2++) {
      double dx = px(at2, frame) - px(at1, frame) ;
      double dy = py(at2, frame) - py(at1, frame) ;
      double dz = pz(at2, frame) - pz(at1, frame) ;

      int hcurr = 0 ;
      int kcurr = 0 ;
      int lcurr = 0 ;
      
      if(usePBC[0] && (pbc[0] || pbc[1] || pbc[2])){
        double da = dx*cellInv(0, 0) + dy*cellInv(0, 1) + dz*cellInv(0, 2) ;
        double db = dx*cellInv(1, 0) + dy*cellInv(1, 1) + dz*cellInv(1, 2) ;
        double dc = dx*cellInv(2, 0) + dy*cellInv(2, 1) + dz*cellInv(2, 2) ;
        if(pbc[0]) {
          // Miller indices of atom 2
          hcurr = floor(da + 0.5) ;
          // Wrap distances
          da = da - hcurr ;
        }
        if(pbc[1]) {
          kcurr = floor(db + 0.5) ;
          db = db - kcurr ;
        }
        if(pbc[2]) {
          lcurr = floor(dc + 0.5) ;
          dc = dc - lcurr ;
        }
        dx = da*cell(0, 0) + db*cell(0, 1) + dc*cell(0, 2) ;
        dy = da*cell(1, 0) + db*cell(1, 1) + dc*cell(1, 2) ;
        dz = da*cell(2, 0) + db*cell(2, 1) + dc*cell(2, 2) ;
      }
      double ds = dx*dx + dy*dy + dz*dz ;
      double dbond = safety[0]*(r[at2] + r[at1]) ;
      if(ds < dbond*dbond) {
        atm1.push_back(at1 + 1) ;
        atm2.push_back(at2 + 1) ;
        h.push_back(hcurr) ;
        k.push_back(kcurr) ;
        l.push_back(lcurr) ;
      }
    }
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = atm1,
    Rcpp::Named("atm2") = atm2,
    Rcpp::Named("h") = h,
    Rcpp::Named("k") = k,
    Rcpp::Named("l") = l) ;
}





// Rcpp::DataFrame guessBondsAtoms2( // Check out VMD: BondSearch.C
//     const Rcpp::S4& x,
//     const Rcpp::Nullable< Rcpp::NumericVector > radius = R_NilValue,
//     const Rcpp::NumericVector& safety = 1.2){
//   std::vector< double > r ;
//   if(radius.isNull()){
//     Rcpp::DataFrame atmprop = x.slot("atoms") ;
//     if(atmprop.containsElementNamed("radius")){
//       r = Rcpp::as< std::vector< double > >(atmprop["radius"]) ;
//     } else {
//       Rcpp::Environment PeriodicTableEnv("package:PeriodicTable") ;
//       Rcpp::Function rcov = PeriodicTableEnv["rcov"] ;
//       r = Rcpp::as< std::vector< double > >(rcov(x)) ;
//     }
//   } else {
//     r = Rcpp::as< std::vector< double > >(radius) ;
//     Rcpp::DataFrame atmprop = x.slot("atoms") ;
//     int rsize = r.size() ;
//     if(rsize != atmprop.nrows()){
//       Rcpp::stop("'radius' must be of length ", atmprop.nrows()) ;
//     }
//   }
//   if(safety.size() != 1)
//     Rcpp::stop("'safety' must be a numeric vector of length 1") ;
//   Rcpp::LogicalVector pbc = x.slot("pbc") ;
//   Rcpp::IntegerVector current = x.slot("current") ;
//   int frame = current[0] - 1 ;
//   Rcpp::NumericMatrix xtraj = x.slot("x") ;
//   Rcpp::NumericMatrix ytraj = x.slot("y") ;
//   Rcpp::NumericMatrix ztraj = x.slot("z") ;
// 
//   Rcpp::Environment AtomsEnv("package:Atoms") ;
//   Rcpp::Function getCell = AtomsEnv["cell"] ;
//   Rcpp::NumericMatrix cell = getCell(x) ;
//   
//   Rcpp::Environment BaseEnv("package:base") ;
//   Rcpp::Function solve = BaseEnv["solve"] ;
//   Rcpp::NumericMatrix cellInv = solve(cell) ;
//   
//   // range?
//   
//   for(int atm = 0; atm != px.nrow(); atm++) {
//     // px
//   }
//   
//   
//   
//   return Rcpp::DataFrame::create() ;
//   
//   // return Rcpp::DataFrame::create(
//   //   Rcpp::Named("atm1") = atm1,
//   //   Rcpp::Named("atm2") = atm2,
//   //   Rcpp::Named("h") = h,
//   //   Rcpp::Named("k") = k,
//   //   Rcpp::Named("l") = l) ;
// }

