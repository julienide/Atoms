#ifndef LAMMPS_READER_H
#define LAMMPS_READER_H

#include "AtomsReader.h"
#include "utils.h"

class LAMMPSReader : public AtomsReader
{
public:
  LAMMPSReader(
    Rcpp::CharacterVector files_,
    Rcpp::IntegerVector selection_,
    Rcpp::IntegerVector first_,
    Rcpp::IntegerVector last_,
    Rcpp::IntegerVector stride_) ;
  
  std::string atomStyle = "" ;
  
  int natoms = 0 ;
  int nbonds = 0 ;
  int nangles = 0 ;
  int ndihedrals = 0 ;
  int nimpropers = 0 ;
  
  int bondTypes = 0 ;
  int angleTypes = 0 ;
  int dihedralTypes = 0 ;
  int improperTypes = 0 ;

  // Attributs
  // atomic:	atom-ID atom-type x y z
  // charge:	atom-ID atom-type q x y z
  // full:	atom-ID molecule-ID atom-type q x y z
  // std::vector<int> atmnumb ;
  std::vector<int> resnumb ;
  std::vector<int> atmtype ;
  std::vector<double> charge ;
  
  // Methods
  void setAtomStyle() ;
  void setAtomicProperties() ;
  void readFrame() ;
} ;

#endif // LAMMPS_READER_H
