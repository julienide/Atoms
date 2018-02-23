#include "LAMMPSReader.h"
#include "utils.h"

LAMMPSReader::LAMMPSReader(
  Rcpp::CharacterVector files_,
  Rcpp::IntegerVector selection_,
  Rcpp::IntegerVector first_,
  Rcpp::IntegerVector last_,
  Rcpp::IntegerVector stride_) :
  AtomsReader(files_, selection_, first_, last_, stride_)
{
  
}

void LAMMPSReader::setAtomicProperties()
{
  if(atomStyle == "atomic")
  {
    atoms = Rcpp::DataFrame::create(
      Rcpp::Named("atmtype") = atmtype,
      Rcpp::Named("stringsAsFactors") = false ) ;
  }
  else if(atomStyle == "charge")
  {
    atoms = Rcpp::DataFrame::create(
      Rcpp::Named("atmtype") = atmtype,
      Rcpp::Named("charge") = charge,
      Rcpp::Named("stringsAsFactors") = false ) ;
  }
  else if(atomStyle == "full")
  {
    atoms = Rcpp::DataFrame::create(
      Rcpp::Named("atmtype") = atmtype,
      Rcpp::Named("resnumb") = resnumb,
      Rcpp::Named("charge") = charge,
      Rcpp::Named("stringsAsFactors") = false ) ;
  }
}

void LAMMPSReader::setAtomStyle(){
  int ncol = 0 ;
  std::istringstream lineReader ;
  lineReader.str(line) ;
  std::string col ;
  while(!lineReader.fail())
  {
    lineReader >> col ;
    if(lineReader.fail() || col.substr(0, 1) == "#")
    {
      break ;
    }
    ncol++ ;
  }
  lineReader.clear() ;
  if(ncol == 5){
    atomStyle = "atomic";
  }
  else if(ncol == 6)
  {
    atomStyle = "charge";
  }
  else if(ncol == 7)
  {
    atomStyle = "full";
  }
  else
  {
    Rcpp::stop("Unrecognized atom style. It should be one of: atomic, charge or full") ;
  }
}

void LAMMPSReader::readFrame(){
  frame++ ;
  if(skipFrame())
  {
    // A LAMMPS file can contain only one single frame
    // Then directly skip file
    fileReader.setstate(std::ios::eofbit) ;
  }
  else
  {
    double xlo = 0.0 ; double xhi = 0.0 ;
    double ylo = 0.0 ; double yhi = 0.0 ;
    double zlo = 0.0 ; double zhi = 0.0 ;
    double xy = 0.0 ; double xz = 0.0 ; double yz = 0.0 ;
    std::istringstream lineReader ;
    nextLine() ; // Title line
    while(!fileReader.eof())
    {
      nextLine() ;
      if(line.find("atoms") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> natoms ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      }
      else if(line.find("bonds") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> nbonds ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      }
      else if(line.find("angles") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> nangles ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      }
      else if(line.find("dihedrals") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> ndihedrals ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      }
      else if(line.find("impropers") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> nimpropers ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("bond" ) != std::string::npos &&
                line.find("types") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> bondTypes ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("angle") != std::string::npos &&
                line.find("types") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> angleTypes ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("dihedral") != std::string::npos &&
                line.find("types") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> dihedralTypes ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("improper") != std::string::npos &&
                line.find("types") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> improperTypes ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("xlo") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> xlo >> xhi ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("ylo") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> ylo >> yhi ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("zlo") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> zlo >> zhi ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("xy") != std::string::npos)
      {
        lineReader.str(line) ;
        lineReader >> xy >> xz >> yz ;
        if(lineReader.fail()){
          errUnexpectedFormat() ;
        }
        lineReader.clear() ;
      } else if(line.find("Atoms") != std::string::npos)
      {
        if(natoms == 0)
        {
          Rcpp::stop("Unspecified number of atoms") ;
        }
        else
        {
          while(atom != natoms)
          {
            nextLine() ;
            line = trim(line) ;
            if(!line.empty()){
              atom++ ;
              if(isFirstAtom())
              {
                nframe++ ;
              }
              if(allAtoms() || (selAtom != selection.end() && atom == *selAtom))
              {
                if(atom == 1){
                  setAtomStyle() ;
                }
                lineReader.str(line) ;
                double anumb ;
                double atype ;
                double x_, y_, z_ ;
                if(atomStyle == "atomic")
                {
                  lineReader >> anumb >> atype >> x_ >> y_ >> z_ ;
                  if(lineReader.fail()){
                    errUnexpectedFormat() ;
                  }
                }
                else if(atomStyle == "charge")
                {
                  double q ;
                  lineReader >> anumb >> atype >> q >> x_ >> y_ >> z_ ;
                  if(lineReader.fail()){
                    errUnexpectedFormat() ;
                  }
                  charge.push_back(q) ;
                }
                else if(atomStyle == "full")
                {
                  double q ;
                  double rnumb ;
                  lineReader >> anumb >> rnumb >> atype >> q >> x_ >> y_ >> z_ ;
                  if(lineReader.fail()){
                    errUnexpectedFormat() ;
                  }
                  charge.push_back(q) ;
                  resnumb.push_back(rnumb) ;
                }
                lineReader.clear() ;
                atmtype.push_back(atype) ;
                x.push_back(x_) ;
                y.push_back(y_) ;
                z.push_back(z_) ;
              }
              if(!allAtoms())
              {
                if(!isFirstFrame() && selAtom == selection.end())
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
        }
      }
      else if(line.find("Bonds") != std::string::npos)
      {
        if(isFirstFrame())
        {
          int B = 0 ;
          while(B != nbonds && !fileReader.eof())
          {
            nextLine() ;
            line = trim(line) ;
            if(!line.empty())
            {
              B++ ;
              int seq, type, atm1, atm2 ;
              lineReader.str(line) ;
              lineReader >> seq >> type >> atm1 >> atm2 ;
              if(lineReader.fail()){
                errUnexpectedFormat() ;
              }
              lineReader.clear() ;
              Batm1.push_back(atm1) ;
              Batm2.push_back(atm2) ;
            }
          }
          if(B != nbonds)
          {
            Rcpp::stop("Missing Bonds") ;
          }
        }
      }
      else if(line.find("Angles") != std::string::npos)
      {
        if(isFirstFrame())
        {
          int A = 0 ;
          while(A != nangles && !fileReader.eof())
          {
            nextLine() ;
            line = trim(line) ;
            if(!line.empty())
            {
              A++ ;
              int seq, type, atm1, atm2, atm3 ;
              lineReader.str(line) ;
              lineReader >> seq >> type >> atm1 >> atm2 >> atm3 ;
              if(lineReader.fail()){
                errUnexpectedFormat() ;
              }
              lineReader.clear() ;
              Aatm1.push_back(atm1) ;
              Aatm2.push_back(atm2) ;
              Aatm3.push_back(atm3) ;
            }
          }
          if(A != nangles)
          {
            Rcpp::stop("Missing Angles") ;
          }
        }
      }
      else if(line.find("Dihedrals") != std::string::npos)
      {
        if(isFirstFrame())
        {
          int D = 0 ;
          while(D != ndihedrals && !fileReader.eof())
          {
            nextLine() ;
            line = trim(line) ;
            if(!line.empty())
            {
              D++ ;
              int seq, type, atm1, atm2, atm3, atm4 ;
              lineReader.str(line) ;
              lineReader >> seq >> type >> atm1 >> atm2 >> atm3 >> atm4 ;
              if(lineReader.fail()){
                errUnexpectedFormat() ;
              }
              lineReader.clear() ;
              Datm1.push_back(atm1) ;
              Datm2.push_back(atm2) ;
              Datm3.push_back(atm3) ;
              Datm4.push_back(atm4) ;
            }
          }
          if(D != ndihedrals)
          {
            Rcpp::stop("Missing Dihedrals") ;
          }
        }
      }
      else if(line.find("Impropers") != std::string::npos)
      {
        if(isFirstFrame())
        {
          int I = 0 ;
          while(I != nimpropers && !fileReader.eof())
          {
            nextLine() ;
            line = trim(line) ;
            if(!line.empty())
            {
              I++ ;
              int seq, type, atm1, atm2, atm3, atm4 ;
              lineReader.str(line) ;
              lineReader >> seq >> type >> atm1 >> atm2 >> atm3 >> atm4 ;
              if(lineReader.fail()){
                errUnexpectedFormat() ;
              }
              lineReader.clear() ;
              Iatm1.push_back(atm1) ;
              Iatm2.push_back(atm2) ;
              Iatm3.push_back(atm3) ;
              Iatm4.push_back(atm4) ;
            }
          }
          if(I != nimpropers)
          {
            Rcpp::stop("Missing Impropers") ;
          }
        }
      }
      else if(line.find("Bond") != std::string::npos &&
              line.find("Coeffs") != std::string::npos)
      {
        
      }
      else if(line.find("Angle") != std::string::npos &&
              line.find("Coeffs") != std::string::npos)
      {
        
      }
      else if(line.find("Dihedral") != std::string::npos &&
              line.find("Coeffs") != std::string::npos)
      {
        
      }
      else if(line.find("Improper") != std::string::npos &&
              line.find("Coeffs") != std::string::npos)
      {
        
      }
    }
    if(!allAtoms() && (selAtom != selection.end())){
      errBadSelection(*selAtom) ;
    }
    if(xhi - xlo <= 0.0){
      Rcpp::stop("xhi - xlo <= 0") ;
    }
    if(yhi - ylo <= 0.0){
      Rcpp::stop("yhi - ylo <= 0") ;
    }
    if(zhi - zlo <= 0.0){
      Rcpp::stop("zhi - zlo <= 0") ;
    }
    a.push_back(xhi - xlo) ; a.push_back(0.0      ) ; a.push_back(0.0      ) ;
    b.push_back(xy       ) ; b.push_back(yhi - ylo) ; b.push_back(0.0      ) ;
    c.push_back(xz       ) ; c.push_back(yz       ) ; c.push_back(zhi - zlo) ;
    pbc[0] = true ; pbc[1] = true ; pbc[2] = true ;
  }
}
