#include <Rcpp.h>
#include "guessDihedrals.h"

//' @rdname guessTopology
//' @export
// [[Rcpp::export(name = "guessDihedrals.Atoms")]]
Rcpp::DataFrame guessDihedralsAtoms(const Rcpp::S4& x){
  Rcpp::DataFrame angles = x.slot("angles") ;
  Rcpp::IntegerVector atm1 = angles["atm1"] ;
  Rcpp::IntegerVector atm2 = angles["atm2"] ;
  Rcpp::IntegerVector atm3 = angles["atm3"] ;

  std::vector<int> Datm1 ;
  std::vector<int> Datm2 ;
  std::vector<int> Datm3 ;
  std::vector<int> Datm4 ;

  for(int a1 = 0; a1 != atm1.size(); ++a1) {
    for(int a2 = a1 + 1; a2 < atm1.size(); ++a2) {
      if(atm1[a1] == atm2[a2] && atm2[a1] == atm3[a2]){
        bool rev = (atm1[a1] > atm2[a1]) ;
        rev = rev || (atm1[a1] == atm2[a1] && atm1[a2] > atm3[a1]) ;
        if(rev){
          Datm1.push_back(atm3[a1]) ;
          Datm2.push_back(atm2[a1]) ;
          Datm3.push_back(atm1[a1]) ;
          Datm4.push_back(atm1[a2]) ;
        } else {
          Datm1.push_back(atm1[a2]) ;
          Datm2.push_back(atm1[a1]) ;
          Datm3.push_back(atm2[a1]) ;
          Datm4.push_back(atm3[a1]) ;
        }
        // if(atm1[a1] < atm2[a1]){
        //   Datm1.push_back(atm1[a2]) ;
        //   Datm2.push_back(atm1[a1]) ;
        //   Datm3.push_back(atm2[a1]) ;
        //   Datm4.push_back(atm3[a1]) ;
        // } else {
        //   Datm1.push_back(atm3[a1]) ;
        //   Datm2.push_back(atm2[a1]) ;
        //   Datm3.push_back(atm1[a1]) ;
        //   Datm4.push_back(atm1[a2]) ;
        // }
      } else if(atm1[a1] == atm2[a2] && atm2[a1] == atm1[a2]){
        bool rev = (atm1[a1] > atm2[a1]) ;
        rev = rev || (atm1[a1] == atm2[a1] && atm3[a2] > atm3[a1]) ;
        if(rev){
          Datm1.push_back(atm3[a1]) ;
          Datm2.push_back(atm2[a1]) ;
          Datm3.push_back(atm1[a1]) ;
          Datm4.push_back(atm3[a2]) ;
        } else {
          Datm1.push_back(atm3[a2]) ;
          Datm2.push_back(atm1[a1]) ;
          Datm3.push_back(atm2[a1]) ;
          Datm4.push_back(atm3[a1]) ;
        }
        // if(atm1[a1] < atm2[a1]){
        //   Datm1.push_back(atm3[a2]) ;
        //   Datm2.push_back(atm1[a1]) ;
        //   Datm3.push_back(atm2[a1]) ;
        //   Datm4.push_back(atm3[a1]) ;
        // } else {
        //   Datm1.push_back(atm3[a1]) ;
        //   Datm2.push_back(atm2[a1]) ;
        //   Datm3.push_back(atm1[a1]) ;
        //   Datm4.push_back(atm3[a2]) ;
        // }
      } else if(atm2[a1] == atm1[a2] && atm3[a1] == atm2[a2]){
        bool rev = (atm2[a1] > atm3[a1]) ;
        rev = rev || (atm2[a1] == atm3[a1] && atm1[a1] > atm3[a2]) ;
        if(rev){
          Datm1.push_back(atm3[a2]) ;
          Datm2.push_back(atm3[a1]) ;
          Datm3.push_back(atm2[a1]) ;
          Datm4.push_back(atm1[a1]) ;
        } else {
            Datm1.push_back(atm1[a1]) ;
            Datm2.push_back(atm2[a1]) ;
            Datm3.push_back(atm3[a1]) ;
            Datm4.push_back(atm3[a2]) ;
        }
        // if(atm2[a1] < atm3[a1]){
        //   Datm1.push_back(atm1[a1]) ;
        //   Datm2.push_back(atm2[a1]) ;
        //   Datm3.push_back(atm3[a1]) ;
        //   Datm4.push_back(atm3[a2]) ;
        // } else {
        //   Datm1.push_back(atm3[a2]) ;
        //   Datm2.push_back(atm3[a1]) ;
        //   Datm3.push_back(atm2[a1]) ;
        //   Datm4.push_back(atm1[a1]) ;
        // }
      } else if(atm2[a1] == atm3[a2] && atm3[a1] == atm2[a2]){
        bool rev = (atm2[a1] > atm3[a1]) ;
        rev = rev || (atm2[a1] == atm3[a1] && atm1[a1] > atm1[a2]) ;
        if(rev){
          Datm1.push_back(atm1[a2]) ;
          Datm2.push_back(atm3[a1]) ;
          Datm3.push_back(atm2[a1]) ;
          Datm4.push_back(atm1[a1]) ;
        } else {
          Datm1.push_back(atm1[a1]) ;
          Datm2.push_back(atm2[a1]) ;
          Datm3.push_back(atm3[a1]) ;
          Datm4.push_back(atm1[a2]) ;
        }
        // if(atm2[a1] < atm3[a1]){
        //   Datm1.push_back(atm1[a1]) ;
        //   Datm2.push_back(atm2[a1]) ;
        //   Datm3.push_back(atm3[a1]) ;
        //   Datm4.push_back(atm1[a2]) ;
        // } else {
        //   Datm1.push_back(atm1[a2]) ;
        //   Datm2.push_back(atm3[a1]) ;
        //   Datm3.push_back(atm2[a1]) ;
        //   Datm4.push_back(atm1[a1]) ;
        // }
      }
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Datm1,
    Rcpp::Named("atm2") = Datm2,
    Rcpp::Named("atm3") = Datm3,
    Rcpp::Named("atm4") = Datm4) ;
}
