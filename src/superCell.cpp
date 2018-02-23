// #include "superCell.h"
// 
// //' @rdname superCell
// //' @export
// //' @noRd
// // [[Rcpp::export(name = ".superCellAtoms")]]
// Rcpp::S4 superCellAtoms(
//     const Rcpp::S4& x,
//     const Rcpp::IntegerVector& aInds,
//     const Rcpp::IntegerVector& bInds,
//     const Rcpp::IntegerVector& cInds)
// {
//   Rcpp::DataFrame atoms = x.slot("atoms") ;
//   
//   Rcpp::LogicalVector pbc = x.slot("pbc") ;
//   Rcpp::NumericMatrix a = x.slot("a") ;
//   Rcpp::NumericMatrix b = x.slot("b") ;
//   Rcpp::NumericMatrix c = x.slot("c") ;
// 
//   Rcpp::Environment AtomsEnv("package:Atoms") ;
//   Rcpp::Function fractional = AtomsEnv["fractional"] ;
//   Rcpp::DataFrame frac = fractional(x) ;
//   Rcpp::NumericVector xFrac = frac["x"] ;
//   Rcpp::NumericVector yFrac = frac["y"] ;
//   Rcpp::NumericVector zFrac = frac["z"] ;
//   
//   if(!pbc[0] && aInds.size() != 1){
//     Rcpp::stop("No PBC along a-axis") ;
//   }
//   if(!pbc[1] && bInds.size() != 1){
//     Rcpp::stop("No PBC along b-axis") ;
//   }
//   if(!pbc[2] && cInds.size() != 1){
//     Rcpp::stop("No PBC along c-axis") ;
//   }
//   
//   int natom = frac.nrows() ;
//   int ncell = aInds.size()*bInds.size()*cInds.size() ;
//   Rcpp::NumericMatrix xout(natom*ncell, 1) ;
//   Rcpp::NumericMatrix yout(natom*ncell, 1) ;
//   Rcpp::NumericMatrix zout(natom*ncell, 1) ;
//   int row = 0 ;
//   for(Rcpp::IntegerVector::const_iterator ai = aInds.begin(); ai != aInds.end(); ai++)
//   {
//     for(Rcpp::IntegerVector::const_iterator bi = bInds.begin(); bi != bInds.end(); bi++)
//     {
//       for(Rcpp::IntegerVector::const_iterator ci = cInds.begin(); ci != cInds.end(); ci++)
//       {
//         Rcpp::NumericVector::iterator xFracIt = xFrac.begin() ;
//         Rcpp::NumericVector::iterator yFracIt = yFrac.begin() ;
//         Rcpp::NumericVector::iterator zFracIt = zFrac.begin() ;
//         for(int atm = 0; atm < natom; atm++){
//           double aTrans = (*xFracIt + *ai) ;
//           double bTrans = (*yFracIt + *bi) ;
//           double cTrans = (*zFracIt + *ci) ;
//           xout(row, 0) = aTrans*a(0, 0) + bTrans*b(0, 0) + cTrans*c(0, 0) ;
//           yout(row, 0) = aTrans*a(1, 0) + bTrans*b(1, 0) + cTrans*c(1, 0) ;
//           zout(row, 0) = aTrans*a(2, 0) + bTrans*b(2, 0) + cTrans*c(2, 0) ;
//           xFracIt++ ;
//           yFracIt++ ;
//           zFracIt++ ;
//           row++ ;
//         }
//       }
//     }
//   }
//   
//   Rcpp::S4 Obj("Atoms") ;
//   Obj.slot("pbc") = pbc ;
//   Obj.slot("a") = a*aInds.size() ;
//   Obj.slot("b") = b*bInds.size() ;
//   Obj.slot("c") = c*cInds.size() ;
//   Obj.slot("current") = Rcpp::IntegerVector::create(1) ;
//   Obj.slot("x") = xout ;
//   Obj.slot("y") = yout ;
//   Obj.slot("z") = zout ;
//   
//   return Obj ;
// }