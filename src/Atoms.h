#ifndef ATOMS_H
#define ATOMS_H

#include <Rcpp.h>

class Atoms
{
public:
  Atoms() ;
  
  // Attributs
  Rcpp::DataFrame atoms ;
  std::vector< int > Batm1, Batm2 ;
  std::vector< int > Aatm1, Aatm2, Aatm3 ;
  std::vector< int > Datm1, Datm2, Datm3, Datm4 ;
  std::vector< int > Iatm1, Iatm2, Iatm3, Iatm4 ;
  int natom = 0 ; int nframe = 0 ;
  std::vector< bool > pbc = {false, false, false} ;
  std::vector< double > a, b, c ;
  std::vector< double > x, y, z ;
  std::vector< double > vx, vy, vz ;
  std::vector< double > fx, fy, fz ;
  
  // Methods
  virtual void setAtomicProperties() =0 ;
  Rcpp::S4 asS4() ;
} ;

#endif // ATOMS_H
