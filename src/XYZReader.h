#ifndef XYZ_READER_H
#define XYZ_READER_H

#include "AtomsReader.h"

class XYZReader : public AtomsReader
{
public:
  XYZReader(
    Rcpp::CharacterVector files_,
    Rcpp::IntegerVector selection_,
    Rcpp::IntegerVector first_,
    Rcpp::IntegerVector last_,
    Rcpp::IntegerVector stride_) ;
  
  // Attributs
  std::vector< std::string > atmtype ;
  
  // Methods
  void setAtomicProperties() ;
  void readFrame() ;
} ;

#endif // XYZ_READER_H
