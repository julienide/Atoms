#include "XYZWriter.h"

//' @export
//' @rdname writeAtoms
//' @noRd
// [[Rcpp::export(name = ".XYZWriter")]]
void XYZWriter(
    Rcpp::S4 x,
    Rcpp::CharacterVector file){
  
  bool multi = file.size() != 1 ;
  
  Rcpp::DataFrame atoms = x.slot("atoms") ;
  Rcpp::CharacterVector atmtype ;
  if(atoms.containsElementNamed("atmtype") &&
     TYPEOF(atoms["atmtype"]) == STRSXP){
    atmtype = atoms["atmtype"] ;
  } else {
    atmtype.push_back("Xx") ;
  }
  Rcpp::Environment PeriodicTable("package:PeriodicTable") ;
  Rcpp::Function getSymb = PeriodicTable["symb"] ;
  Rcpp::CharacterVector symb = getSymb(atmtype) ;
  
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  unsigned int natom = px.nrow() ;
  unsigned int nframe = px.ncol() ;

  for(unsigned int frame = 0; frame < nframe; frame++){
    std::string filename ;
    std::ofstream fileWriter ;
    if(multi){
      filename = Rcpp::as<std::string>(file[frame]) ;
    } else {
      filename = Rcpp::as<std::string>(file[0]) ;
    }
    if(frame == 0 || multi){
      fileWriter.open(filename, std::fstream::out) ;
    } else {
      fileWriter.open(filename, std::fstream::out | std::fstream::app) ;
    }
    fileWriter << natom << std::endl << std::endl ;
    for(unsigned int atom = 0; atom < natom; atom++){
      std::ostringstream lineWriter ;
      lineWriter << std::left ;
      lineWriter << std::setw(3) << symb[atom] << " " ;
      lineWriter << std::right ;
      lineWriter << std::fixed << std::setw(12) << std::setprecision(6) << px(atom, frame) << " " ;
      lineWriter << std::fixed << std::setw(12) << std::setprecision(6) << py(atom, frame) << " " ;
      lineWriter << std::fixed << std::setw(12) << std::setprecision(6) << pz(atom, frame) << " " ;
      fileWriter << lineWriter.str() << std::endl ;
    }
  }
}