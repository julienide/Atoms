#include <Rcpp.h>
#include "guessAngles.h"

//' @rdname guessTopology
//' @export
// [[Rcpp::export(name = "guessAngles.Atoms")]]
Rcpp::DataFrame guessAnglesAtoms(const Rcpp::S4& x){
  Rcpp::DataFrame bonds = x.slot("bonds") ;
  Rcpp::IntegerVector atm1 = bonds["atm1"] ;
  Rcpp::IntegerVector atm2 = bonds["atm2"] ;
  
  std::vector<int> Aatm1 ;
  std::vector<int> Aatm2 ;
  std::vector<int> Aatm3 ;
  
  for(int b1 = 0; b1 != atm1.size(); ++b1) {
    for(int b2 = b1 + 1; b2 < atm1.size(); ++b2) {
      if(atm1[b1] == atm1[b2]){
        Aatm2.push_back(atm1[b1]) ;
        // if(atm2[b1] < atm2[b2]){
          Aatm1.push_back(atm2[b1]) ;
          Aatm3.push_back(atm2[b2]) ;
          // c(atm2[b1], atm1[b1], atm2[b2])
        // } else {
          // Aatm1.push_back(atm2[b1]) ;
          // Aatm3.push_back(atm2[b2]) ;
          // c(atm2[b2], atm1[b1], atm2[b1])
        // }
      } else if(atm1[b1] == atm2[b2]){
        Aatm2.push_back(atm1[b1]) ;
        // if(atm2[b1] < atm1[b2]){
          Aatm1.push_back(atm2[b1]) ;
          Aatm3.push_back(atm1[b2]) ;
          // c(atm2[b1], atm1[b1], atm1[b2])
        // } else {
          // Aatm1.push_back(atm1[b2]) ;
          // Aatm3.push_back(atm2[b1]) ;
          // c(atm1[b2], atm1[b1], atm2[b1])
        // }
      } else if(atm2[b1] == atm1[b2]){
        Aatm2.push_back(atm2[b1]) ;
        // if(atm1[b1] < atm2[b2]){
          Aatm1.push_back(atm1[b1]) ;
          Aatm3.push_back(atm2[b2]) ;
          // c(atm1[b1], atm2[b1], atm2[b2])
        // } else {
          // Aatm1.push_back(atm2[b2]) ;
          // Aatm3.push_back(atm1[b1]) ;
          // c(atm2[b2], atm2[b1], atm1[b1])
        // }
      } else if(atm2[b1] == atm2[b2]){
        Aatm2.push_back(atm2[b1]) ;
        // if(atm1[b1] < atm1[b2]){
          Aatm1.push_back(atm1[b1]) ;
          Aatm3.push_back(atm1[b2]) ;
          // c(atm1[b1], atm2[b1], atm1[b2])
        // } else {
          // Aatm1.push_back(atm1[b2]) ;
          // Aatm3.push_back(atm1[b1]) ;
          // c(atm1[b2], atm2[b1], atm1[b1])
        // }
      }
    }
  }
  
  // Rcpp::DataFrame Obj = Rcpp::DataFrame::create(
  //   Rcpp::Named("atm1") = Aatm1,
  //   Rcpp::Named("atm2") = Aatm2,
  //   Rcpp::Named("atm3") = Aatm3) ;
  // 
  // Rcpp::Environment BaseEnv("package:base") ;
  // Rcpp::Function order = BaseEnv["order"] ;
  // Rcpp::IntegerVector O = order(Obj) ;
  // 
  // for(Rcpp::IntgerVector::Iterator a = O.begin(); a != O.end(); ++O){
  //   
  // }
  
  
  
  // return sort(Obj) ;
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("atm1") = Aatm1,
    Rcpp::Named("atm2") = Aatm2,
    Rcpp::Named("atm3") = Aatm3) ;
}
