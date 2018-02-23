#include "RESReader.h"
#include "utils.h"

RESReader::RESReader(
  Rcpp::CharacterVector files_,
  Rcpp::IntegerVector selection_,
  Rcpp::IntegerVector first_,
  Rcpp::IntegerVector last_,
  Rcpp::IntegerVector stride_) :
  AtomsReader(files_, selection_, first_, last_, stride_)
{
  
}

void RESReader::setAtomicProperties()
{
  atoms = Rcpp::DataFrame::create(
    Rcpp::Named("atmname") = atmname,
    Rcpp::Named("atmtype") = atmtype,
    Rcpp::Named("stringsAsFactors") = false ) ;
}

void RESReader::readFrame(){
  frame++ ;
  if(skipFrame())
  {
    // A RES file can contain only one single frame
    // Then directly skip file
    fileReader.setstate(std::ios::eofbit) ;
  }
  else
  {
    bool CELL = false ;
    bool SFAC = false ;
    std::vector< std::string > sfac ;
    while(!fileReader.eof())
    {
      if(line.length())
      {
        std::istringstream lineReader ; lineReader.str(line) ;
        std::string record ; lineReader >> record ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        if(record == "ABIN"){
          
        } else if(record == "ACTA"){
          
        } else if(record == "AFIX"){
          
        } else if(record == "ANIS"){
          
        } else if(record == "ANSC"){
          
        } else if(record == "ANSR"){
          
        } else if(record == "BASF"){
          
        } else if(record == "BIND"){
          
        } else if(record == "BLOC"){
          
        } else if(record == "BOND"){
          
        } else if(record == "BUMP"){
          
        } else if(record == "CELL"){
          if(CELL)
          {
            errUnexpectedRecord("CELL") ;
          }
          double lambda ;
          lineReader >> lambda >> aLength >> bLength >> cLength >> cosA >> cosB >> cosG ;
          if(lineReader.fail()){
            errUnexpectedFormat() ;
          }
          cosA = cos(cosA*pi/180.0) ; cosB = cos(cosB*pi/180.0) ; cosG = cos(cosG*pi/180.0) ;
          double sinG = sqrt(1 - pow(cosG, 2.0)) ;
          a.push_back(aLength) ;
          a.push_back(0.0) ;
          a.push_back(0.0) ;
          b.push_back(bLength*cosG) ;
          b.push_back(bLength*sinG) ;
          b.push_back(0.0) ;
          c.push_back(cLength*cosB) ;
          c.push_back(cLength*(cosA - cosB*cosG)/sinG) ;
          c.push_back(cLength*sqrt(1 - pow(cosA, 2.0) - pow(cosB , 2.0) - pow(cosG, 2.0) + 2*(cosA*cosB*cosG))/sinG) ;
          pbc[0] = true ; pbc[1] = true ; pbc[2] = true ;
          CELL = true ;
        } else if(record == "CGLS"){
          
        } else if(record == "CHIV"){
          
        } else if(record == "CONF"){
          
        } else if(record == "CONN"){
          
        } else if(record == "DAMP"){
          
        } else if(record == "DANG"){
          
        } else if(record == "DEFS"){
          
        } else if(record == "DELU"){
          
        } else if(record == "DFIX"){
          
        } else if(record == "DISP"){
          
        } else if(record == "EADP"){
          
        } else if(record == "END"){
          
        } else if(record == "EQIV"){
          
        } else if(record == "EXTI"){
          
        } else if(record == "EXYZ"){
          
        } else if(record == "FEND"){
          
        } else if(record == "FLAT"){
          
        } else if(record == "FMAP"){
          
        } else if(record == "FRAG"){
          
        } else if(record == "FREE"){
          
        } else if(record == "FVAR"){
          
        } else if(record == "GRID"){
          
        } else if(record == "HFIX"){
          
        } else if(record == "HKLF"){
          
        } else if(record == "HTAB"){
          
        } else if(record == "ISOR"){
          
        } else if(record == "LATT"){
          
        } else if(record == "LAUE"){
          
        } else if(record == "LIST"){
          
        } else if(record == "L.S."){
          
        } else if(record == "MERG"){
          
        } else if(record == "MORE"){
          
        } else if(record == "MOVE"){
          
        } else if(record == "MPLA"){
          
        } else if(record == "NCSY"){
          
        } else if(record == "NEUT"){
          
        } else if(record == "OMIT"){
          
        } else if(record == "PART"){
          
        } else if(record == "PLAN"){
          
        } else if(record == "PRIG"){
          
        } else if(record == "REM"){
          
        } else if(record == "RESI"){
          
        } else if(record == "RIGU"){
          
        } else if(record == "RTAB"){
          
        } else if(record == "SADI"){
          
        } else if(record == "SAME"){
          
        } else if(record == "SFAC"){
          if(SFAC) // SFAC has already been read
          {
            errUnexpectedRecord("SFAC") ;
          }
          if(!sfac.size())
          {
            std::string atype ; 
            lineReader >> atype ;
            while(!lineReader.fail())
            {
              sfac.push_back(atype) ;
              lineReader >> atype ;
            }
          }
        } else if(record == "SHEL"){
          
        } else if(record == "SIMU"){
          
        } else if(record == "SIZE"){
          
        } else if(record == "SPEC"){
          
        } else if(record == "STIR"){
          
        } else if(record == "SUMP"){
          
        } else if(record == "SWAT"){
          
        } else if(record == "SYMM"){
          
        } else if(record == "TEMP"){
          
        } else if(record == "TITL"){
          
        } else if(record == "TWIN"){
          
        } else if(record == "TWST"){
          
        } else if(record == "UNIT"){
          
        } else if(record == "WGHT"){
          
        } else if(record == "WIGL"){
          
        } else if(record == "WPDB"){
          
        } else if(record == "XNPD"){
          
        } else if(record == "ZERR"){
          
        } else {
          atom++ ;
          if(isFirstAtom())
          {
            nframe++ ;
          }
          if(allAtoms() || (selAtom != selection.end() && atom == *selAtom))
          {
            if(!CELL)
            {
              Rcpp::stop("'CELL' was not specified") ;
            }
            if(!sfac.size())
            {
              Rcpp::stop("'SFAC' was not specified or is empty") ;
            }
            int anumb ; double xp, yp, zp ;
            lineReader >> anumb >> xp >> yp >> zp ;
            if(lineReader.fail())
            {
              errUnexpectedFormat() ;
            }
            if(static_cast<int>(sfac.size()) < anumb)
            {
              Rcpp::stop("Invalid scattering factor number: ", anumb) ;
            }
            if(isFirstFrame())
            {
              atmname.push_back(record) ;
              atmtype.push_back(sfac[anumb - 1]) ;
            }
            x.push_back(a[0]*xp + b[0]*yp + c[0]*zp) ;
            y.push_back(a[1]*xp + b[1]*yp + c[1]*zp) ;
            z.push_back(a[2]*xp + b[2]*yp + c[2]*zp) ;
            if(!allAtoms())
            {
              if(selAtom == selection.end())
              {
                // To exit the file and move to the next file/frame
                fileReader.setstate(std::ios::eofbit) ;
                break ;
              }
              else
              {
                selAtom++ ;
              }
            }
          }
        }
        lineReader.clear() ;
      }
      nextLine() ;
    }
    if(!allAtoms() && (selAtom != selection.end())){
      errBadSelection(*selAtom) ;
    }
  }
}
