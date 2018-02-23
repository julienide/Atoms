#include "RESWriter.h"

const double pi = std::atan(1.0)*4 ;

//' @export
//' @rdname writeAtoms
//' @noRd
// [[Rcpp::export(name = ".RESWriter")]]
void RESWriter(
    Rcpp::S4 x,
    Rcpp::CharacterVector file){
  
  Rcpp::DataFrame atoms = x.slot("atoms") ;
  std::vector< std::string > atmname ;
  if(atoms.containsElementNamed("atmname") &&
     TYPEOF(atoms["atmname"]) == STRSXP){
    atmname = Rcpp::as<std::vector< std::string >>(atoms["atmname"]) ;
  } else {
    atmname.push_back("Xx") ;
  }
  std::vector< std::string > atmtype ;
  if(atoms.containsElementNamed("atmtype") &&
     TYPEOF(atoms["atmtype"]) == STRSXP){
    atmtype = Rcpp::as<std::vector< std::string >>(atoms["atmtype"]) ;
  } else {
    atmtype.push_back("Xx") ;
  }
  
  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  unsigned int natom = px.nrow() ;
  unsigned int nframe = px.ncol() ;
  
  for(unsigned int frame = 0; frame < nframe; frame++){
    std::string filename ;
    std::ofstream fileWriter ;
    filename = Rcpp::as<std::string>(file[frame]) ;
    fileWriter.open(filename, std::fstream::out) ;
    fileWriter << "TITLE " << filename << std::endl ;
    double aLength = sqrt(
      a(0, frame)*a(0, frame) +
      a(1, frame)*a(1, frame) +
      a(2, frame)*a(2, frame)) ;
    double bLength = sqrt(
      b(0, frame)*b(0, frame) +
      b(1, frame)*b(1, frame) +
      b(2, frame)*b(2, frame)) ;
    double cLength = sqrt(
      c(0, frame)*c(0, frame) +
      c(1, frame)*c(1, frame) +
      c(2, frame)*c(2, frame)) ;
    double alpha = 
      b(0, frame)*c(0, frame) +
      b(1, frame)*c(1, frame) +
      b(2, frame)*c(2, frame) ;
    double beta  = 
      c(0, frame)*a(0, frame) +
      c(1, frame)*a(1, frame) +
      c(2, frame)*a(2, frame) ;
    double gamma = 
      a(0, frame)*b(0, frame) +
      a(1, frame)*b(1, frame) +
      a(2, frame)*b(2, frame) ;
    alpha = acos(alpha/(bLength*cLength))*180.0/pi ;
    beta  = acos(beta /(cLength*aLength))*180.0/pi ;
    gamma = acos(gamma/(aLength*bLength))*180.0/pi ;
    std::ostringstream lineWriter ;
    lineWriter << "CELL "  ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << aLength ;
    lineWriter << std::setw(1) << " " ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << bLength ;
    lineWriter << std::setw(1) << " " ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << cLength ;
    lineWriter << std::setw(1) << " " ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << alpha ;
    lineWriter << std::setw(1) << " " ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << beta ;
    lineWriter << std::setw(1) << " " ;
    lineWriter << std::fixed << std::setw(8) << std::setprecision(4) << gamma ;
    fileWriter << lineWriter.str() << std::endl ;
    fileWriter << "LATT -1" << std::endl ;
    lineWriter.str("") ;
    lineWriter << "SFAC " ;
    std::vector< std::string > sfac(atmtype.begin(), atmtype.end()) ;
    std::sort(sfac.begin(), sfac.end()) ;
    std::vector< std::string >::iterator it = std::unique(sfac.begin(), sfac.end()) ;
    std::map<std::string, unsigned int> sfacMap ;
    for(std::vector< std::string >::iterator sf = sfac.begin(); sf != it; sf++){
      sfacMap[*sf] = std::distance(sfac.begin(), sf) + 1 ;
      lineWriter << " " << *sf ;
    }
    fileWriter << lineWriter.str() << std::endl ;
    for(unsigned int atom = 0; atom < natom; atom++){
      lineWriter.str("") ;
      lineWriter << std::left << std::setw(5) << atmname[atom].substr(0, 5) << " " ;
      lineWriter << std::right << std::setw(3) << sfacMap[atmtype[atom]] << " " ;
      lineWriter << std::fixed << std::setw(10) << std::setprecision(5) << px(atom, frame) << " " ;
      lineWriter << std::fixed << std::setw(10) << std::setprecision(5) << py(atom, frame) << " " ;
      lineWriter << std::fixed << std::setw(10) << std::setprecision(5) << pz(atom, frame) << " " ;
      lineWriter << "   1.00000   0.00000" ;
      fileWriter << lineWriter.str() << std::endl ;
    }
    fileWriter << "END" << std::endl ;
  }
}