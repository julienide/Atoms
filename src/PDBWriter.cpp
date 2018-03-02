#include "PDBWriter.h"

const double pi = std::atan(1.0)*4 ;

//' @export
//' @rdname writeAtoms
//' @noRd
// [[Rcpp::export(name = ".PDBWriter")]]
void PDBWriter(
  Rcpp::S4 x,
  Rcpp::CharacterVector file){

  bool multi = file.size() != 1 ;

  Rcpp::DataFrame atoms = x.slot("atoms") ;
  std::vector< std::string > recname ;
  if(atoms.containsElementNamed("recname") &&
     TYPEOF(atoms["recname"]) == STRSXP){
    recname = Rcpp::as< std::vector< std::string > >(atoms["recname"]) ;
  } else {
    recname.push_back("ATOM") ;
  }
  std::vector< std::string > atmname ;
  if(atoms.containsElementNamed("atmname") &&
     TYPEOF(atoms["atmname"]) == STRSXP){
    atmname = Rcpp::as< std::vector< std::string > >(atoms["atmname"]) ;
  } else {
    atmname.push_back("Xx") ;
  }
  std::vector< std::string > alt ;
  if(atoms.containsElementNamed("alt") &&
     TYPEOF(atoms["alt"]) == STRSXP){
    alt = Rcpp::as< std::vector< std::string > >(atoms["alt"]) ;
  } else {
    alt.push_back(" ") ;
  }
  std::vector< std::string > resname ;
  if(atoms.containsElementNamed("resname") &&
     TYPEOF(atoms["resname"]) == STRSXP){
    resname = Rcpp::as< std::vector< std::string > >(atoms["resname"]) ;
  } else {
    resname.push_back("UNK") ;
  }
  std::vector< std::string > chain ;
  if(atoms.containsElementNamed("chain") &&
     TYPEOF(atoms["chain"]) == STRSXP){
    chain = Rcpp::as< std::vector< std::string > >(atoms["chain"]) ;
  } else {
    chain.push_back(" ") ;
  }
  std::vector< int > resnumb ;
  if(atoms.containsElementNamed("resnumb") &&
     TYPEOF(atoms["resnumb"]) == INTSXP){
    resnumb = Rcpp::as< std::vector< int > >(atoms["resnumb"]) ;
  } else {
    resnumb.push_back(1) ;
  }
  std::vector< std::string > insert ;
  if(atoms.containsElementNamed("insert") &&
     TYPEOF(atoms["insert"]) == STRSXP){
    insert = Rcpp::as< std::vector< std::string > >(atoms["insert"]) ;
  } else {
    insert.push_back(" ") ;
  }
  std::vector< double > occ ;
  if(atoms.containsElementNamed("occ") &&
     TYPEOF(atoms["occ"]) == REALSXP){
    occ = Rcpp::as< std::vector< double > >(atoms["occ"]) ;
  } else {
    occ.push_back(0.0) ;
  }
  std::vector< double > temp ;
  if(atoms.containsElementNamed("temp") &&
     TYPEOF(atoms["temp"]) == REALSXP){
    temp = Rcpp::as< std::vector< double > >(atoms["temp"]) ;
  } else {
    temp.push_back(0.0) ;
  }
  std::vector< std::string > segname ;
  if(atoms.containsElementNamed("segname") &&
     TYPEOF(atoms["segname"]) == STRSXP){
    segname = Rcpp::as< std::vector< std::string > >(atoms["segname"]) ;
  } else {
    segname.push_back("") ;
  }
  std::vector< std::string > atmtype ;
  if(atoms.containsElementNamed("atmtype") &&
     TYPEOF(atoms["atmtype"]) == STRSXP){
    atmtype = Rcpp::as< std::vector< std::string > >(atoms["atmtype"]) ;
  } else {
    atmtype.push_back("Xx") ;
  }
  Rcpp::Environment PeriodicTable("package:PeriodicTable") ;
  Rcpp::Function getSymb = PeriodicTable["symb"] ;
  std::vector< std::string > symb = Rcpp::as< std::vector< std::string > >(getSymb(atmtype)) ;

  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  unsigned int natom = px.nrow() ;
  unsigned int nframe = px.ncol() ;
  Rcpp::DataFrame bonds = x.slot("bonds") ;

  // std::string filename = Rcpp::as<std::string>(file[0]) ;
  // std::ofstream fileWriter(filename) ;

  Progress pb(nframe, true) ;
  for(unsigned int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    std::string filename ;
    std::ofstream fileWriter ;
    if(multi){
      filename = Rcpp::as<std::string>(file[frame]) ;
    } else {
      filename = Rcpp::as<std::string>(file[0]) ;
    }
    if(frame == 0 || multi){
      fileWriter.open(filename, std::fstream::out) ;
    } else {
      fileWriter.open(filename, std::fstream::out | std::fstream::app) ;
    }
    // MODEL
    if(nframe > 1 && !multi){
      std::ostringstream lineWriter ;
      lineWriter << std::setw(6) << std::left << "MODEL" ;
      lineWriter << "    " ;
      lineWriter << std::setw(4) << std::right << frame + 1 ;
      fileWriter << lineWriter.str() << "\n" ;
    }
    // CRYST1
    if(pbc[0] && pbc[1] && pbc[2]){
      double aLength = sqrt(
        a(0, frame)*a(0, frame) +
        a(1, frame)*a(1, frame) +
        a(2, frame)*a(2, frame)) ;
      double bLength = sqrt(
        b(0, frame)*b(0, frame) +
        b(1, frame)*b(1, frame) +
        b(2, frame)*b(2, frame)) ;
      double cLength = sqrt(
        c(0, frame)*c(0, frame) +
        c(1, frame)*c(1, frame) +
        c(2, frame)*c(2, frame)) ;
      double alpha =
        b(0, frame)*c(0, frame) +
        b(1, frame)*c(1, frame) +
        b(2, frame)*c(2, frame) ;
      double beta  =
        c(0, frame)*a(0, frame) +
        c(1, frame)*a(1, frame) +
        c(2, frame)*a(2, frame) ;
      double gamma =
        a(0, frame)*b(0, frame) +
        a(1, frame)*b(1, frame) +
        a(2, frame)*b(2, frame) ;
      alpha = acos(alpha/(bLength*cLength))*180.0/pi ;
      beta  = acos(beta /(cLength*aLength))*180.0/pi ;
      gamma = acos(gamma/(aLength*bLength))*180.0/pi ;
      std::ostringstream lineWriter ;
      lineWriter << "CRYST1" ;
      lineWriter << std::fixed << std::setw(9) << std::setprecision(3) << aLength ;
      lineWriter << std::fixed << std::setw(9) << std::setprecision(3) << bLength ;
      lineWriter << std::fixed << std::setw(9) << std::setprecision(3) << cLength ;
      lineWriter << std::fixed << std::setw(7) << std::setprecision(2) << alpha ;
      lineWriter << std::fixed << std::setw(7) << std::setprecision(2) << beta ;
      lineWriter << std::fixed << std::setw(7) << std::setprecision(2) << gamma ;
      lineWriter << std::setw(16) << " P 1           1" ;
      fileWriter << lineWriter.str() << std::endl ;
    }
    for(unsigned int atom = 0; atom < natom; atom++){
      // ATOM & HETATM
      //
      // COLUMNS        DATA  TYPE    FIELD        DEFINITION
      // -------------------------------------------------------------------------------------
      //  1 -  6        Record name   "ATOM  "
      //  7 - 11        Integer       serial       Atom  serial number.
      // 13 - 16        Atom          name         Atom name.
      // 17             Character     altLoc       Alternate location indicator.
      // 18 - 20        Residue name  resName      Residue name.
      // 22             Character     chainID      Chain identifier.
      // 23 - 26        Integer       resSeq       Residue sequence number.
      // 27             AChar         iCode        Code for insertion of residues.
      // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
      // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
      // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
      // 55 - 60        Real(6.2)     occupancy    Occupancy.
      // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
      // 77 - 78        LString(2)    element      Element symbol, right-justified.
      // 79 - 80        LString(2)    charge       Charge  on the atom.
      std::ostringstream lineWriter ;
      if(recname.size() == natom){
        lineWriter << std::left << std::setw(6) << recname[atom].substr(0, 6) ;
      } else {
        lineWriter << std::left << std::setw(6) << recname[0].substr(0, 6) ;
      }
      if(atom + 1 <= 99999){
        lineWriter << std::right << std::setw(5) << atom + 1 ;
      } else if(atom + 1 <= 1048575){
        lineWriter << std::right << std::setw(5) << std::hex << atom + 1 ;
      } else {
        lineWriter << std::setw(5) << "*****" ;
      }
      if(atom + 1 > 99999){
      } else {
      }
      lineWriter << std::setw(1) << " " ;
      if(atmname.size() == natom){
        // if(trim(atmname[atom]).length() == 3){
        //   " XXX"
        // }
        // if(trim(atmname[atom]).length() == 4){
        //   "XXXX"
        // }
        lineWriter << std::left << std::setw(4) << atmname[atom].substr(0, 4) ;
      } else {
        lineWriter << std::left << std::setw(4) << atmname[0].substr(0, 4) ;
      }
      if(alt.size() == natom){
        lineWriter << std::setw(1) << alt[atom].substr(0, 1) ;
      } else {
        lineWriter << std::setw(1) << alt[0].substr(0, 1) ;
      }
      if(resname.size() == natom){
        lineWriter << std::setw(4) << resname[atom].substr(0, 4) ;
      } else {
        lineWriter << std::setw(4) << resname[0].substr(0, 4) ;
      }
      if(chain.size() == natom){
        lineWriter << std::setw(1) << chain[atom].substr(0, 1) ;
      } else {
        lineWriter << std::setw(1) << chain[0].substr(0, 1) ;
      }
      if(resnumb.size() == natom){
        lineWriter << std::fixed << std::right << std::setw(4) << resnumb[atom] ;
      } else {
        lineWriter << std::fixed << std::right << std::setw(4) << resnumb[0] ;
      }
      if(insert.size() == natom){
        lineWriter << std::setw(1) << insert[atom].substr(0, 1) ;
      } else {
        lineWriter << std::setw(1) << insert[0].substr(0, 1) ;
      }
      lineWriter << "   " ;
      lineWriter << std::fixed << std::setw(8) << std::setprecision(3) << px(atom, frame) ;
      lineWriter << std::fixed << std::setw(8) << std::setprecision(3) << py(atom, frame) ;
      lineWriter << std::fixed << std::setw(8) << std::setprecision(3) << pz(atom, frame) ;
      if(occ.size() == natom){
        lineWriter << std::fixed << std::setw(6) << std::setprecision(2) << occ[atom] ;
      } else {
        lineWriter << std::fixed << std::setw(6) << std::setprecision(2) << occ[0] ;
      }
      if(temp.size() == natom){
        lineWriter << std::fixed << std::setw(6) << std::setprecision(2) << temp[atom] ;
      } else {
        lineWriter << std::fixed << std::setw(6) << std::setprecision(2) << temp[0] ;
      }
      lineWriter << std::setw(6) << "      " ;
      if(segname.size() == natom){
        lineWriter << std::left << std::setw(4) << segname[atom].substr(0, 4) ;
      } else {
        lineWriter << std::left << std::setw(4) << segname[0].substr(0, 4) ;
      }
      if(symb.size() == natom){
        lineWriter << std::right << std::setw(2) << symb[atom].substr(0, 2) ;
      } else {
        lineWriter << std::right << std::setw(2) << symb[0].substr(0, 2) ;
      }
      fileWriter << lineWriter.str() << std::endl ;
    }
    // ENDMDL
    if(nframe > 1 && !multi){
      std::ostringstream lineWriter ;
      lineWriter << std::setw(6) << std::left << "ENDMDL" ;
      fileWriter << lineWriter.str() << std::endl ;
    }
    if(frame + 1 == nframe || multi){
      // CONECT
      if(bonds.nrows() != 0){
        Rcpp::IntegerVector atm1 = bonds["atm1"] ;
        Rcpp::IntegerVector atm2 = bonds["atm2"] ;
        int nbond = atm1.size() ;
        int curAtm1 = -1 ;
        int natm2 = 0 ;
        int bond = 0 ;
        while(bond < nbond){
          curAtm1 = atm1[bond] ;
          natm2 = 0 ;
          std::ostringstream lineWriter ;
          lineWriter << "CONECT" ;
          lineWriter << std::right << std::setw(5) << atm1[bond] ;
          while(atm1[bond] == curAtm1){
            lineWriter << std::right << std::setw(5) << atm2[bond] ;
            bond++ ;
            natm2++ ;
            if(natm2 == 4){
              break ;
            }
          }
          fileWriter << lineWriter.str() << std::endl ;
        }
      }
      // END
      fileWriter << "END" << "\n" ;
    }
    pb.increment() ;
  }
}
