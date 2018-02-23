#include "rdf.h"

//' @rdname rdf
//' @export
// [[Rcpp::export(name = ".rdfAtoms")]]
Rcpp::DataFrame rdfAtoms(
    const Rcpp::S4& x,
    const Rcpp::IntegerVector& sel1,
    const Rcpp::IntegerVector& sel2,
    const Rcpp::IntegerVector& resnumb,
    const Rcpp::NumericVector& cutoff,
    const Rcpp::NumericVector& interval,
    const std::vector< bool >& type){
  double _cutoff = cutoff[0] ;
  double _interval = interval[0] ;
  
  int dim = type[0] + type[1] + type[2] ;
  double Vfac = 1.0 ;
  if(dim > 1)
  {
    Vfac = Vfac*std::atan(1.0)*4 ; // multiply by pi
  }
  if(dim == 3)
  {
    Vfac = Vfac*4.0/3.0 ; // 4/3 for the volume of a sphere
  }
  
  bool lowertri = true ;
  Rcpp::IntegerVector::const_iterator S1 = sel1.begin() ;
  Rcpp::IntegerVector::const_iterator S2 = sel2.begin() ;
  for(S1 = sel1.begin() ; S1 != sel1.end() ; S1++)
  {
    if(S2 == sel2.end() || S1 != S2)
    {
      lowertri = false ;
      break ;
    }
    S2++ ;
  }

  Rcpp::LogicalVector pbc = x.slot("pbc") ;
  Rcpp::NumericMatrix a = x.slot("a") ;
  Rcpp::NumericMatrix b = x.slot("b") ;
  Rcpp::NumericMatrix c = x.slot("c") ;
  Rcpp::NumericMatrix px = x.slot("x") ;
  Rcpp::NumericMatrix py = x.slot("y") ;
  Rcpp::NumericMatrix pz = x.slot("z") ;
  int nframe = a.ncol() ;
  
  std::vector< int > hkl(3, 0) ;
  for(int i = 0; i < 3; i++){
    if(pbc[i]){
      hkl[i] = 1 ;
    }
  }
  
  int nbin = (int) _cutoff/_interval ;
  std::vector< double > r(nbin) ;
  std::vector< double > gInter(nbin) ;
  std::vector< double > gIntra(nbin) ;
  std::vector< double > cnInter(nbin) ;
  std::vector< double > cnIntra(nbin) ;

  Rcpp::Environment baseEnv("package:base") ;
  Rcpp::Function solve = baseEnv["solve"] ;

  Progress pb(nframe, true) ;
  for(int frame = 0; frame < nframe; frame++){
    if(Progress::check_abort()){
      Rcpp::stop("Interuption") ;
    }
    
    // Rcpp::NumericVector curPx = px(Rcpp::_, frame) ;
    // Rcpp::NumericVector curPy = px(Rcpp::_, frame) ;
    // Rcpp::NumericVector curPz = px(Rcpp::_, frame) ;
    
    Rcpp::NumericMatrix cell(3, 3) ;
    cell(0, 0) = a(0, frame) ; cell(0, 1) = b(0, frame) ; cell(0, 2) = c(0, frame) ;
    cell(1, 0) = a(1, frame) ; cell(1, 1) = b(1, frame) ; cell(1, 2) = c(1, frame) ;
    cell(2, 0) = a(2, frame) ; cell(2, 1) = b(2, frame) ; cell(2, 2) = c(2, frame) ;
    Rcpp::NumericMatrix cellInv = solve(cell) ;
    
    unsigned int ntotal = 0 ;
    std::vector< double > inter(nbin) ;
    std::vector< double > intra(nbin) ;
    for(Rcpp::IntegerVector::const_iterator S1 = sel1.begin(); S1 != sel1.end(); S1++){
      int A1 = *S1 - 1 ;
      for(Rcpp::IntegerVector::const_iterator S2 = sel2.begin(); S2 != sel2.end(); S2++){
        if(lowertri && (*S2 >= *S1)){
          break ;
        }
        int A2 = *S2 - 1 ;
        std::vector< double > d(3, 0.0) ;
        if(type[0])
        {
          d[0] = px(A2, frame) - px(A1, frame) ;
        }
        if(type[1])
        {
          d[1] = py(A2, frame) - py(A1, frame) ;
        }
        if(type[2])
        {
          d[2] = pz(A2, frame) - pz(A1, frame) ;
        }
        if(pbc[0] || pbc[1] || pbc[2])
        {
          applyPBC2(d[0], d[1], d[2], cell, cellInv, pbc) ; // Takes some time...
        }
        double dn = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) ;
        if(dn <= _cutoff){
          dn = dn/_interval ;
          int bin = (int) dn ;
          if(resnumb[A1] == resnumb[A2]){
            intra[bin] = intra[bin] + 1.0 ;
          } else {
            inter[bin] = inter[bin] + 1.0 ;
          }
        }

        // // RDF with neighboring cells:
        // // Much longer but allow for calculating the RDF up to cutoff < L instead of L/2
        // for(int h = -hkl[0]; h <= hkl[0]; h++){
        //   for(int k = -hkl[1]; k <= hkl[1]; k++){
        //     for(int l = -hkl[2]; l <= hkl[2]; l++){
        //       double dn =
        //         (d[0] + h*a(0, frame) + k*b(0, frame) + l*c(0, frame))*
        //         (d[0] + h*a(0, frame) + k*b(0, frame) + l*c(0, frame)) +
        //         (d[1] + h*a(1, frame) + k*b(1, frame) + l*c(1, frame))*
        //         (d[1] + h*a(1, frame) + k*b(1, frame) + l*c(1, frame)) +
        //         (d[2] + h*a(2, frame) + k*b(2, frame) + l*c(2, frame))*
        //         (d[2] + h*a(2, frame) + k*b(2, frame) + l*c(2, frame)) ;
        //       dn = sqrt(dn) ;
        //       if(dn <= _cutoff){
        //         dn = dn/_interval ;
        //         int bin = (int) dn ;
        //         if(resnumb[A1] == resnumb[A2]){
        //           intra[bin] = intra[bin] + 1.0 ;
        //         } else {
        //           inter[bin] = inter[bin] + 1.0 ;
        //         }
        //       }
        //     }
        //   }
        // }

        ntotal++ ;
      }
    }
    // Cell Volume/Area/Length: V = |a.(b x c)| ; A = |a x b| ; L = |a|
    std::vector< double > u(3, 0.0), v(3, 0.0), w(3, 0.0) ;
    u[0] = 1.0 ; v[1] = 1.0 ; w[2] = 1.0 ;
    if(type[0])
    {
      u[0] = a(0, frame) ;
      u[1] = a(1, frame) ;
      u[2] = a(2, frame) ;
    }
    if(type[1])
    {
      v[0] = b(0, frame) ;
      v[1] = b(1, frame) ;
      v[2] = b(2, frame) ;
    }
    if(type[2])
    {
      w[0] = c(0, frame) ;
      w[1] = c(1, frame) ;
      w[2] = c(2, frame) ;
    }
    double V =
      u[0] * (v[1]*w[2] - v[2]*w[1]) +
      u[1] * (v[2]*w[0] - v[0]*w[2]) +
      u[2] * (v[0]*w[1] - v[1]*w[0]) ;
    
    for(int bin = 0; bin < nbin; bin++)
    {
      if(ntotal != 0)
      {
        inter[bin] = inter[bin]/ntotal ;
        intra[bin] = intra[bin]/ntotal ;
      }
      inter[bin] = inter[bin]/nframe ;
      intra[bin] = intra[bin]/nframe ;
      gInter[bin] = gInter[bin] + inter[bin]*V ;
      gIntra[bin] = gIntra[bin] + intra[bin]*V ;
      cnInter[bin] = cnInter[bin] + inter[bin] ;
      cnIntra[bin] = cnIntra[bin] + intra[bin] ;
    }
    pb.increment() ;
  }
  
  double prevInter = 0.0 ;
  double prevIntra = 0.0 ;
  for(int bin = 0 ; bin < nbin ; bin++)
  {
    // Distance
    r[bin] = bin*_interval + _interval/2.0 ;

    // Normalization by bin's volume
    double R1 = (r[bin] - _interval/2L) ;
    double R2 = (r[bin] + _interval/2L) ;
    double binVol = Vfac*(pow(R2, dim) - pow(R1, dim)) ;
    gInter[bin] = gInter[bin]/binVol ;
    gIntra[bin] = gIntra[bin]/binVol ;

    // Cumulative sum to calculate coordination numbers of sel1.
    cnInter[bin] = (prevInter + cnInter[bin]*sel2.size()) ;
    cnIntra[bin] = (prevIntra + cnIntra[bin]*sel2.size()) ;
    prevInter = cnInter[bin] ;
    prevIntra = cnIntra[bin] ;
  }

  Rcpp::DataFrame Obj = Rcpp::DataFrame::create(
    Rcpp::Named("r") = r,
    Rcpp::Named("gInter") = gInter,
    Rcpp::Named("gIntra") = gIntra,
    Rcpp::Named("cnInter") = cnInter,
    Rcpp::Named("cnIntra") = cnIntra) ;
  return Obj ;
}
