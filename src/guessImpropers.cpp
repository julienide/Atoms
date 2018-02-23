#include <Rcpp.h>
#include "guessImpropers.h"

// Works by assuming that if (1,2,3) is an angle, and 2 & 4 are bonded,
// then (2, 1, 3, 4) must be an improper dihedral.
// ie the improper dihedral is the angle between the planes formed by
// (1, 2, 3) and (1, 3, 4)

//' @rdname guessTopology
//' @export
// [[Rcpp::export(name = "guessImpropers.Atoms")]]
Rcpp::DataFrame guessImpropersAtoms(const Rcpp::S4& x){
  Rcpp::DataFrame angles = x.slot("angles") ;
  Rcpp::IntegerVector atm1 = angles["atm1"] ;
  Rcpp::IntegerVector atm2 = angles["atm2"] ;
  Rcpp::IntegerVector atm3 = angles["atm3"] ;
  
  std::vector<int> Iatm1 ;
  std::vector<int> Iatm2 ;
  std::vector<int> Iatm3 ;
  std::vector<int> Iatm4 ;
  
  for(int a1 = 0; a1 != atm1.size(); ++a1) {
    for(int a2 = a1 + 1; a2 < atm1.size(); ++a2) {
      if(atm2[a1] == atm2[a2]){
        if(atm1[a1] == atm1[a2]){
          if(atm1[a1] < atm3[a1]){
            Iatm1.push_back(atm2[a1]) ;
            Iatm2.push_back(atm1[a1]) ;
            Iatm3.push_back(atm3[a1]) ;
            Iatm4.push_back(atm3[a2]) ;
          } else {
            Iatm1.push_back(atm3[a2]) ;
            Iatm2.push_back(atm1[a1]) ;
            Iatm3.push_back(atm3[a1]) ;
            Iatm4.push_back(atm2[a1]) ;
          }
        } else if(atm1[a1] == atm3[a2]){
          if(atm1[a1] < atm3[a1]){
            Iatm1.push_back(atm2[a1]) ;
            Iatm2.push_back(atm1[a1]) ;
            Iatm3.push_back(atm3[a1]) ;
            Iatm4.push_back(atm1[a2]) ;
          } else {
            Iatm1.push_back(atm1[a2]) ;
            Iatm2.push_back(atm1[a1]) ;
            Iatm3.push_back(atm3[a1]) ;
            Iatm4.push_back(atm2[a1]) ;
          }
        } else if(atm3[a1] == atm1[a2]){
          if(atm3[a1] < atm1[a1]){
            Iatm1.push_back(atm2[a1]) ;
            Iatm2.push_back(atm3[a1]) ;
            Iatm3.push_back(atm1[a1]) ;
            Iatm4.push_back(atm3[a2]) ;
          } else {
            Iatm1.push_back(atm3[a2]) ;
            Iatm2.push_back(atm3[a1]) ;
            Iatm3.push_back(atm1[a1]) ;
            Iatm4.push_back(atm2[a1]) ;
          }
        } else if(atm3[a1] == atm3[a2]){
          if(atm3[a1] < atm1[a1]){
            Iatm1.push_back(atm2[a1]) ;
            Iatm2.push_back(atm3[a1]) ;
            Iatm3.push_back(atm1[a1]) ;
            Iatm4.push_back(atm1[a2]) ;
          } else {
            Iatm1.push_back(atm1[a2]) ;
            Iatm2.push_back(atm3[a1]) ;
            Iatm3.push_back(atm1[a1]) ;
            Iatm4.push_back(atm2[a1]) ;
          }
        }
      }
    }
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Iatm1,
    Rcpp::Named("atm2") = Iatm2,
    Rcpp::Named("atm3") = Iatm3,
    Rcpp::Named("atm4") = Iatm4) ;
}
