#ifndef PDB_WRITER_H
#define PDB_WRITER_H

#include <Rcpp.h>
#include <fstream>
#include <progress.hpp>

void PDBWriter(
  Rcpp::S4 x,
  Rcpp::CharacterVector file) ;

#endif // PDB_WRITER_H
