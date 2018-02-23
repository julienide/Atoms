#ifndef RES_WRITER_H
#define RES_WRITER_H

#include <Rcpp.h>
#include <fstream>

void RESWriter(
  Rcpp::S4 x,
  Rcpp::CharacterVector file) ;

#endif // RES_WRITER_H
