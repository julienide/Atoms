#include "AtomsReader.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

AtomsReader::AtomsReader(
  Rcpp::CharacterVector files_,
  Rcpp::IntegerVector selection_,
  Rcpp::IntegerVector first_,
  Rcpp::IntegerVector last_,
  Rcpp::IntegerVector stride_)
{
  files = Rcpp::as< std::vector< std::string > >(files_) ;
  selection = Rcpp::as< std::vector< int > >(selection_) ;
  first = first_[0] ; last = last_[0] ; stride = stride_[0] ;
}

std::string AtomsReader::getFilename()
{
  std::string filename = *file ;
  return filename ;
}

bool AtomsReader::multiFile()
{
  return (files.size() > 1) ;
}

bool AtomsReader::allAtoms()
{
  return (selection[0] == -1) ;
}

bool AtomsReader::isLastAtom()
{
  return (atom == natom) ;
}

bool AtomsReader::isFirstAtom()
{
  return (atom == 1) ;
}

bool AtomsReader::emptyFrame()
{
  return (atom <= 0) ;
}

bool AtomsReader::allFrames()
{
  return (last == -1) ;
}

// bool AtomsReader::isEmpty()
// {
//   return (nframe == 0) ;
// }

bool AtomsReader::isFirstFrame()
{
  return (nframe == 1) ;
}

bool AtomsReader::isLastFrame()
{
  return (!allFrames()) && (frame == last) ;
}

bool AtomsReader::skipFrame()
{
  return (frame < first) || ((frame - first)%stride != 0) ;
}

void AtomsReader::open()
{
  fileReader.open(getFilename(), std::ios::in) ;
  if(!fileReader)
  {
    errMissingFile() ;
  }
}

void AtomsReader::close()
{
  fileReader.close() ;
}

void AtomsReader::setFirstFrame()
{
  natom = x.size() ;
  F1_natom = natom ;
  F1_pbc[0] = pbc[0] ;
  F1_pbc[1] = pbc[1] ;
  F1_pbc[2] = pbc[2] ;
}

void AtomsReader::checkFrame()
{
  int natmCur = x.size()/nframe ;
  if((natmCur != F1_natom) || (x.size()%nframe != 0))
  {
    errInvalidNumberOfAtoms(nframe) ;
  }
  if(pbc[0] != F1_pbc[0] || pbc[1] != F1_pbc[1] || pbc[2] != F1_pbc[2])
  {
    errInvalidPBC(nframe) ;
  }
}

void AtomsReader::initFrame()
{
  aLength = 1.0 ; bLength = 1.0 ; cLength = 1.0 ;
  cosA = 0.0 ; cosB = 0.0 ; cosG = 0.0 ;
  atom = 0 ; line = "" ; selAtom = selection.begin() ;
}

void AtomsReader::nextLine()
{
  std::getline(fileReader, line) ; lineNumber++ ;
}

void AtomsReader::readFile()
{
  lineNumber = 0 ;
  while(!fileReader.eof() && !isLastFrame())
  {
    initFrame() ;
    readFrame() ;
    if(!emptyFrame()){
      if(isFirstFrame()){
        setFirstFrame() ;
      } else {
        checkFrame() ;
      }
    }
  }
}

void AtomsReader::readAllFiles()
{
  if(multiFile()) // Then read one frame per file
  {
    if(allFrames())
    {
      last = files.size() ;
    }
    else
    {
      files.erase(files.begin()+last, files.end()) ;
    }
  }
  
  Progress pb(files.size(), multiFile()) ;
  for(file = files.begin(); file != files.end(); file++)
  {
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    open() ;
    readFile() ;
    close() ;
    pb.increment() ;
    if(isLastFrame()){
      break ;
    }
  }
  if(!allFrames() && frame != last)
  {
    std::string msg = "Less frames read than expected. Last frame read: " +
      std::to_string(frame) ;
    Rcpp::warning(msg) ;
  }
}




// Error messages
void AtomsReader::errMissingFile()
{
  std::string msg = "File '" + getFilename() + "' is missing" ;
  Rcpp::stop(msg) ;
}

std::string AtomsReader::inFileAtLine()
{
  std::string msg =
    " in file '" + getFilename() + "' " +
    "at line: " + std::to_string(lineNumber) ;
  return msg ;
}

void AtomsReader::errUnexpectedFormat()
{
  std::string msg = "Unexpected format" + inFileAtLine() ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errUnexpectedEndOfSection(std::string section)
{
  std::string msg = "Unexpected end of " + section + " section " + inFileAtLine() ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errUnexpectedEndOfFile()
{
  std::string msg = "Unexpected end of file " + getFilename() ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errUnexpectedRecord(std::string record)
{
  std::string msg = "Unexpected " + record + " record " + inFileAtLine() ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errBadSelection(int index)
{
  std::string msg =
    "Atomic selection contains indexes out of range: " + std::to_string(index) ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errInvalidNumberOfAtoms(int nframe)
{
  std::string msg = "Frame " + std::to_string(nframe) + " has an invalid number of atoms" ;
  Rcpp::stop(msg) ;
}

void AtomsReader::errInvalidPBC(int nframe)
{
  std::string msg = "Frame " + std::to_string(nframe) + " has invalid PBC" ;
  Rcpp::stop(msg) ;
}
