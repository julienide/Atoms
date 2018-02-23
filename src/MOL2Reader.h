#ifndef MOL2_READER_H
#define MOL2_READER_H

#include "AtomsReader.h"

class MOL2Reader : public AtomsReader
{
  public:
    MOL2Reader(
      Rcpp::CharacterVector files_,
      Rcpp::IntegerVector selection_,
      Rcpp::IntegerVector first_,
      Rcpp::IntegerVector last_,
      Rcpp::IntegerVector stride_) ;
  
  // Attributs
  std::vector< int > atmnumb ;
  std::vector< std::string > atmname ;
  std::vector< std::string > atmtype ;
  std::vector< int > resnumb ;
  std::vector< std::string > resname ;
  std::vector< double > charge ;
  
  std::istringstream lineReader ;
  
  // Methods
  void setAtomicProperties() ;
  void readFrame() ;
  bool readATOM() ;
} ;

#endif // MOL2_READER_H
