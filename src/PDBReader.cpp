#include "PDBReader.h"
#include "utils.h"

PDBReader::PDBReader(
  Rcpp::CharacterVector files_,
  Rcpp::IntegerVector selection_,
  Rcpp::IntegerVector first_,
  Rcpp::IntegerVector last_,
  Rcpp::IntegerVector stride_) :
  AtomsReader(files_, selection_, first_, last_, stride_)
{
  
}

void PDBReader::setAtomicProperties()
{
  atoms = Rcpp::DataFrame::create(
    Rcpp::Named("recname") = recname,
    Rcpp::Named("atmnumb") = atmnumb,
    Rcpp::Named("atmname") = atmname,
    Rcpp::Named("alt") = alt,
    Rcpp::Named("resname") = resname,
    Rcpp::Named("chain") = chain,
    Rcpp::Named("resnumb") = resnumb,
    Rcpp::Named("insert") = insert,
    Rcpp::Named("occ") = occ,
    Rcpp::Named("temp") = temp,
    Rcpp::Named("segname") = segname,
    Rcpp::Named("atmtype") = atmtype,
    Rcpp::Named("charge") = charge,
    Rcpp::Named("stringsAsFactors") = false ) ;
}

void PDBReader::readFrame(){
  bool MODEL = false ; bool CRYST1 = false ;

  frame++ ;
  while(!fileReader.eof())
  {
    std::string curRecname = rtrim(line.substr(0, 6)) ;
    if(curRecname == "MODEL")
    {
      if(MODEL)
      {
        errUnexpectedRecord("MODEL") ;
      }
      else
      {
        MODEL = true ;
      }
    }
    else if(curRecname == "ENDMDL")
    {
      if(!MODEL)
      {
        errUnexpectedEndOfSection("ENDMDL") ;
      }
      break ;
    }
    else if(curRecname == "CRYST1")
    {
      if(!skipFrame()){
        if(CRYST1)
        {
          errUnexpectedRecord("CRYST1") ;
        }
        else
        {
          CRYST1 = true ;
          line = trim(line) ;
          if(line.length() < 54){
            errUnexpectedFormat() ;
          }
          else
          {
            try {
              aLength = std::stod(line.substr( 6, 9)) ;
              bLength = std::stod(line.substr(15, 9)) ;
              cLength = std::stod(line.substr(25, 9)) ;
              cosA = cos(std::stod(line.substr(33, 7))*pi/180.0) ;
              cosB = cos(std::stod(line.substr(40, 7))*pi/180.0) ;
              cosG = cos(std::stod(line.substr(47, 7))*pi/180.0) ;
            }
            catch(const std::invalid_argument& ia) {
              errUnexpectedFormat() ;
            }
          }
        }
      }
    }
    else if(curRecname == "ATOM" || curRecname == "HETATM")
    {
      if(!skipFrame()){
        atom++ ;
        if(isFirstAtom()){
          nframe++ ;
        }
        if(allAtoms() || (selAtom != selection.end() && atom == *selAtom)){
          line = trim(line) ;
          if(line.length() < 66){
            errUnexpectedFormat() ;
          } else {
            if(isFirstFrame()){
              try {
                recname.push_back(rtrim(curRecname)) ;
                std::string numb = line.substr( 6, 5) ;
                if(atom > 99999){
                  if(numb == "*****")
                  {
                    atmnumb.push_back(atmnumb.back() + 1) ;
                  }
                  else
                  {
                    atmnumb.push_back(std::stoul(numb, nullptr, 16)) ;
                  }
                } else {
                  atmnumb.push_back(std::stoi(numb)) ;
                }
                atmname.push_back(     trim(line.substr(12, 4))) ;
                alt.push_back(         trim(line.substr(16, 1))) ;
                resname.push_back(     trim(line.substr(17, 4))) ;
                chain.push_back(       trim(line.substr(21, 1))) ;
                resnumb.push_back(std::stoi(line.substr(22, 4))) ;
                insert.push_back(      trim(line.substr(26, 1))) ;
              }
              catch(const std::invalid_argument& ia) {
                errUnexpectedFormat() ;
              }
              try{
                occ.push_back(std::stod(line.substr(54, 6))) ;
              }
              catch(const std::invalid_argument& ia) {
                occ.push_back(0.0) ;
              }
              try{
                temp.push_back(std::stod(line.substr(60, 6))) ;
              }
              catch(const std::invalid_argument& ia) {
                temp.push_back(0.0) ;
              }
              try{
                if(line.length() >= 76){ // Sometimes segname are not provided
                  segname.push_back(trim(line.substr(72, 4))) ;
                } else {
                  segname.push_back("") ;
                }
                if(line.length() >= 78){ // Sometimes symbols are not provided
                  atmtype.push_back(trim(line.substr(76, 2))) ;
                } else {
                  atmtype.push_back("") ;
                }
                // if(line.length() >= 80){ // Sometimes charges are not provided
                //   charge = stodNA(line.substr(79, 1) + line.substr(78, 1))
                //   charge.push_back(charge) ;
                // } else {
                //   charge.push_back(0.0) ;
                // }
                charge.push_back(0.0) ;
              }
              catch(const std::invalid_argument& ia) {
                errUnexpectedFormat() ;
              }
            }
            try{
              x.push_back(std::stod(line.substr(30, 8))) ;
              y.push_back(std::stod(line.substr(38, 8))) ;
              z.push_back(std::stod(line.substr(46, 8))) ;
            }
            catch(const std::invalid_argument& ia) {
              errUnexpectedFormat() ;
            }
          }
          if(!allAtoms()){
            if(selAtom == selection.end() && !MODEL && CRYST1){
              break ;
            } else {
              selAtom++ ;
            }
          }
        }
      }
    }
    else if(curRecname == "CONECT")
    {
      line = trim(line) ;
      if(line.length() < 11){
        errUnexpectedFormat() ;
      } else {
        try{
          if(line.length() >= 16){
            Batm1.push_back(std::stoi(line.substr( 6, 5))) ;
            Batm2.push_back(std::stoi(line.substr(11, 5))) ;
          }
          if(line.length() >= 21){
            Batm1.push_back(Batm1.back()) ;
            Batm2.push_back(std::stoi(line.substr(16, 5))) ;
          }
          if(line.length() >= 26){
            Batm1.push_back(Batm1.back()) ;
            Batm2.push_back(std::stoi(line.substr(21, 5))) ;
          }
          if(line.length() >= 31){
            Batm1.push_back(Batm1.back()) ;
            Batm2.push_back(std::stoi(line.substr(26, 5))) ;
          }
        }
        catch(const std::invalid_argument& ia) {
          errUnexpectedFormat() ;
        }
      }
    }
    else if(curRecname == "END")
    {
      break ;
    }
    nextLine() ;
  }
  
  if(emptyFrame())
  {
    if(!skipFrame())
    {
      frame-- ; // Nothing has been read
    }
  } else {
    if(!allAtoms() && (selAtom != selection.end())){
      errBadSelection(*selAtom) ;
    }
    if(CRYST1)
    {
      double sinG = sqrt(1 - pow(cosG, 2.0)) ;
      a.push_back(aLength) ;
      a.push_back(0.0) ;
      a.push_back(0.0) ;
      b.push_back(bLength*cosG) ;
      b.push_back(bLength*sinG) ;
      b.push_back(0.0) ;
      c.push_back(cLength*cosB) ;
      c.push_back(cLength*(cosA - cosB*cosG)/sinG) ;
      c.push_back(cLength*sqrt(1 - pow(cosA, 2.0) - pow(cosB , 2.0) - pow(cosG, 2.0) + 2*cosA*cosB*cosG)/sinG) ;
      pbc[0] = true ;
      pbc[1] = true ;
      pbc[2] = true ;
    }
    else
    {
      if(pbc[0]){
        // If the first frame has PBC then all the frames must have PBC.
        // If there is no CRYST1 record then the first frame's PBC are replicated.
        a.push_back(a[0]) ; a.push_back(a[1]) ; a.push_back(a[2]) ;
        b.push_back(b[0]) ; b.push_back(b[1]) ; b.push_back(b[2]) ;
        c.push_back(c[0]) ; c.push_back(c[1]) ; c.push_back(c[2]) ;
      } else {
        a.push_back(1.0) ; a.push_back(0.0) ; a.push_back(0.0) ;
        b.push_back(0.0) ; b.push_back(1.0) ; b.push_back(0.0) ;
        c.push_back(0.0) ; c.push_back(0.0) ; c.push_back(1.0) ;
      }
    }
  }
  if(files.size() > 1)
  {
    // For multiple PDB file, only the first frame of each file is read
    fileReader.setstate(std::ios::eofbit) ;
  }
}
