#include "XYZReader.h"
#include "utils.h"

XYZReader::XYZReader(
  Rcpp::CharacterVector files_,
  Rcpp::IntegerVector selection_,
  Rcpp::IntegerVector first_,
  Rcpp::IntegerVector last_,
  Rcpp::IntegerVector stride_) :
  AtomsReader(files_, selection_, first_, last_, stride_)
{
  
}

void XYZReader::setAtomicProperties()
{
  atoms = Rcpp::DataFrame::create(
    Rcpp::Named("atmtype") = atmtype,
    Rcpp::Named("stringsAsFactors") = false ) ;
}

void XYZReader::readFrame(){
  frame++ ;
  std::istringstream lineReader ;
  nextLine() ;
  if(fileReader.eof())
  {
    errUnexpectedEndOfFile() ;
  }
  int natm = 0 ;
  lineReader.str(line) ;
  lineReader >> natm ;
  if(lineReader.fail()){
    errUnexpectedFormat() ;
  }
  lineReader.clear() ;
  nextLine() ; // Comment line
  for(atom = 1; atom <= natm; atom++)
  {
    nextLine() ;
    if(fileReader.eof())
    {
      errUnexpectedEndOfFile() ;
    }
    if(!skipFrame())
    {
      if(isFirstAtom())
      {
        nframe++ ;
      }
      if(allAtoms() || (selAtom != selection.end() && atom == *selAtom))
      {
        lineReader.str(line) ;
        std::string atmtype_ ; double x_, y_, z_ ;
        lineReader >> atmtype_ >> x_ >> y_ >> z_ ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
        if(isFirstFrame())
        {
          atmtype.push_back(atmtype_) ;
        }
        x.push_back(x_) ;
        y.push_back(y_) ;
        z.push_back(z_) ;
        if(!allAtoms()){
          if(selAtom == selection.end() && files.size() > 1){
            break ;
          } else {
          selAtom++ ;
          }
        }
      }
    }
  }
  if(skipFrame())
  {
    atom = 0 ;
  }
  else
  {
    if(!allAtoms() && (selAtom != selection.end())){
      errBadSelection(*selAtom) ;
    }
    a.push_back(1.0) ; a.push_back(0.0) ; a.push_back(0.0) ;
    b.push_back(0.0) ; b.push_back(1.0) ; b.push_back(0.0) ;
    c.push_back(0.0) ; c.push_back(0.0) ; c.push_back(1.0) ;
  }
  // Check if next line is EOF
  fileReader.get() ;
  if(!fileReader.eof())
  {
    fileReader.unget() ;
  }
}
