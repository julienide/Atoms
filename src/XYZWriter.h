#ifndef XYZ_WRITER_H
#define XYZ_WRITER_H

#include <Rcpp.h>
#include <fstream>

void XYZWriter(
  Rcpp::S4 x,
  Rcpp::CharacterVector file) ;

#endif // XYZ_WRITER_H
