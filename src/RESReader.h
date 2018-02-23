#ifndef RES_READER_H
#define RES_READER_H

#include "AtomsReader.h"

class RESReader : public AtomsReader
{
public:
  RESReader(
    Rcpp::CharacterVector files_,
    Rcpp::IntegerVector selection_,
    Rcpp::IntegerVector first_,
    Rcpp::IntegerVector last_,
    Rcpp::IntegerVector stride_) ;
  
  // Attributs
  std::vector< std::string > atmname ;
  std::vector< std::string > atmtype ;
  
  // Methods
  void setAtomicProperties() ;
  void readFrame() ;
} ;

#endif // RES_READER_H
