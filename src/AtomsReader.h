#ifndef ATOMS_READER_H
#define ATOMS_READER_H

#include "Atoms.h"
#include <fstream>

const double pi = std::atan(1.0)*4 ;

class AtomsReader : public Atoms
{
public:
  AtomsReader(
    Rcpp::CharacterVector files_,
    Rcpp::IntegerVector selection_,
    Rcpp::IntegerVector first_,
    Rcpp::IntegerVector last_,
    Rcpp::IntegerVector stride_) ;
  std::vector< std::string > files ;
  std::vector< int > selection ;
  int first ; int last ; int stride ;
  
  std::vector< std::string >::iterator file ; // Name of the current file
  std::fstream fileReader ;
  int lineNumber = 0 ; // Number of the current line
  std::string line = "" ; // Current line
  
  double aLength = 1.0 ;
  double bLength = 1.0 ;
  double cLength = 1.0 ;
  double cosA = 0.0 ;
  double cosB = 0.0 ;
  double cosG = 0.0 ;

  int atom = 0 ; // Current atom number
  std::vector< int >::iterator selAtom ; // Next selected atom to read
  
  // Keep track of the current frame over all the files to be read
  int frame = 0 ;
  
  // Save data on the first frame to check the validity of the next frames
  int F1_natom = 0 ;
  std::vector< bool > F1_pbc = {false, false, false} ;
  
  
  // Methods
  std::string getFilename() ; // Return the current file name as a string
  bool multiFile() ; // Whether to read a single file or a batch of files
  bool allAtoms() ; // Whether to read all atoms or not
  bool isLastAtom() ; // Is the current the last one to read
  bool isFirstAtom() ; // Is the current atom the first one
  bool emptyFrame() ; // Is the crrent frame empty
  bool allFrames() ; // Whether to read all the frames or not
  // bool isEmpty() ; // Whether a frame has already been read or not
  bool isFirstFrame() ; // Check if the current frame is the first
  bool isLastFrame() ; // Check if the current frame is the last to be read
  bool skipFrame() ; // Whether to skip the current frame or not
  virtual void setAtomicProperties() =0 ;
  void open() ;
  void close() ;
  void setFirstFrame() ;
  void checkFrame() ;
  void initFrame() ;
  void nextLine() ;
  virtual void readFrame() =0 ;
  void readFile() ;
  void readAllFiles() ;
  
  // Error messages
  void errMissingFile() ;
  std::string inFileAtLine() ;
  void errUnexpectedFormat() ;
  void errUnexpectedEndOfSection(std::string section) ;
  void errUnexpectedEndOfFile() ;
  void errUnexpectedRecord(std::string record) ;
  void errBadSelection(int index) ;
  void errInvalidNumberOfAtoms(int nframe) ;
  void errInvalidPBC(int nframe) ;

} ;

#endif // ATOMS_READER_H
