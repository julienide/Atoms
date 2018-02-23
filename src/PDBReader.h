#ifndef PDB_READER_H
#define PDB_READER_H

#include "AtomsReader.h"

class PDBReader : public AtomsReader
{
public:
  PDBReader(
    Rcpp::CharacterVector files_,
    Rcpp::IntegerVector selection_,
    Rcpp::IntegerVector first_,
    Rcpp::IntegerVector last_,
    Rcpp::IntegerVector stride_) ;

  // Attributs
  std::vector< std::string > recname ;
  std::vector< int > atmnumb ;
  std::vector< std::string > atmname ;
  std::vector< std::string > alt ;
  std::vector< std::string > resname ;
  std::vector< std::string > chain ;
  std::vector< int > resnumb ;
  std::vector< std::string > insert ;
  std::vector< double > occ ;
  std::vector< double > temp ;
  std::vector< std::string > segname ;
  std::vector< std::string > atmtype ;
  std::vector< double > charge ;

  // Methods
  void setAtomicProperties() ;
  void readFrame() ;
} ;

#endif // PDB_READER_H
