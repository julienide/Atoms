// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// LAMMPSParser
Rcpp::S4 LAMMPSParser(Rcpp::CharacterVector files, Rcpp::IntegerVector selection, Rcpp::IntegerVector first, Rcpp::IntegerVector last, Rcpp::IntegerVector stride);
RcppExport SEXP _Atoms_LAMMPSParser(SEXP filesSEXP, SEXP selectionSEXP, SEXP firstSEXP, SEXP lastSEXP, SEXP strideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type files(filesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type last(lastSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type stride(strideSEXP);
    rcpp_result_gen = Rcpp::wrap(LAMMPSParser(files, selection, first, last, stride));
    return rcpp_result_gen;
END_RCPP
}
// MOL2Parser
Rcpp::S4 MOL2Parser(Rcpp::CharacterVector files, Rcpp::IntegerVector selection, Rcpp::IntegerVector first, Rcpp::IntegerVector last, Rcpp::IntegerVector stride);
RcppExport SEXP _Atoms_MOL2Parser(SEXP filesSEXP, SEXP selectionSEXP, SEXP firstSEXP, SEXP lastSEXP, SEXP strideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type files(filesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type last(lastSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type stride(strideSEXP);
    rcpp_result_gen = Rcpp::wrap(MOL2Parser(files, selection, first, last, stride));
    return rcpp_result_gen;
END_RCPP
}
// PDBParser
Rcpp::S4 PDBParser(Rcpp::CharacterVector files, Rcpp::IntegerVector selection, Rcpp::IntegerVector first, Rcpp::IntegerVector last, Rcpp::IntegerVector stride);
RcppExport SEXP _Atoms_PDBParser(SEXP filesSEXP, SEXP selectionSEXP, SEXP firstSEXP, SEXP lastSEXP, SEXP strideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type files(filesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type last(lastSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type stride(strideSEXP);
    rcpp_result_gen = Rcpp::wrap(PDBParser(files, selection, first, last, stride));
    return rcpp_result_gen;
END_RCPP
}
// PDBWriter
void PDBWriter(Rcpp::S4 x, Rcpp::CharacterVector file);
RcppExport SEXP _Atoms_PDBWriter(SEXP xSEXP, SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type file(fileSEXP);
    PDBWriter(x, file);
    return R_NilValue;
END_RCPP
}
// RESParser
Rcpp::S4 RESParser(Rcpp::CharacterVector files, Rcpp::IntegerVector selection, Rcpp::IntegerVector first, Rcpp::IntegerVector last, Rcpp::IntegerVector stride);
RcppExport SEXP _Atoms_RESParser(SEXP filesSEXP, SEXP selectionSEXP, SEXP firstSEXP, SEXP lastSEXP, SEXP strideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type files(filesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type last(lastSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type stride(strideSEXP);
    rcpp_result_gen = Rcpp::wrap(RESParser(files, selection, first, last, stride));
    return rcpp_result_gen;
END_RCPP
}
// RESWriter
void RESWriter(Rcpp::S4 x, Rcpp::CharacterVector file);
RcppExport SEXP _Atoms_RESWriter(SEXP xSEXP, SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type file(fileSEXP);
    RESWriter(x, file);
    return R_NilValue;
END_RCPP
}
// XYZParser
Rcpp::S4 XYZParser(Rcpp::CharacterVector files, Rcpp::IntegerVector selection, Rcpp::IntegerVector first, Rcpp::IntegerVector last, Rcpp::IntegerVector stride);
RcppExport SEXP _Atoms_XYZParser(SEXP filesSEXP, SEXP selectionSEXP, SEXP firstSEXP, SEXP lastSEXP, SEXP strideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type files(filesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type last(lastSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type stride(strideSEXP);
    rcpp_result_gen = Rcpp::wrap(XYZParser(files, selection, first, last, stride));
    return rcpp_result_gen;
END_RCPP
}
// XYZWriter
void XYZWriter(Rcpp::S4 x, Rcpp::CharacterVector file);
RcppExport SEXP _Atoms_XYZWriter(SEXP xSEXP, SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type file(fileSEXP);
    XYZWriter(x, file);
    return R_NilValue;
END_RCPP
}
// comAtoms
Rcpp::S4 comAtoms(const Rcpp::S4& x, const Rcpp::NumericVector& mass, const Rcpp::IntegerVector& factor);
RcppExport SEXP _Atoms_comAtoms(SEXP xSEXP, SEXP massSEXP, SEXP factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass(massSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type factor(factorSEXP);
    rcpp_result_gen = Rcpp::wrap(comAtoms(x, mass, factor));
    return rcpp_result_gen;
END_RCPP
}
// freeVolumeAtoms
Rcpp::NumericVector freeVolumeAtoms(const Rcpp::S4& x, const Rcpp::Nullable< Rcpp::NumericVector > radius);
RcppExport SEXP _Atoms_freeVolumeAtoms(SEXP xSEXP, SEXP radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable< Rcpp::NumericVector > >::type radius(radiusSEXP);
    rcpp_result_gen = Rcpp::wrap(freeVolumeAtoms(x, radius));
    return rcpp_result_gen;
END_RCPP
}
// freeVolume2Atoms
Rcpp::NumericVector freeVolume2Atoms(const Rcpp::S4& x, const Rcpp::Nullable< Rcpp::NumericVector > radius);
RcppExport SEXP _Atoms_freeVolume2Atoms(SEXP xSEXP, SEXP radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable< Rcpp::NumericVector > >::type radius(radiusSEXP);
    rcpp_result_gen = Rcpp::wrap(freeVolume2Atoms(x, radius));
    return rcpp_result_gen;
END_RCPP
}
// guessAnglesAtoms
Rcpp::DataFrame guessAnglesAtoms(const Rcpp::S4& x);
RcppExport SEXP _Atoms_guessAnglesAtoms(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(guessAnglesAtoms(x));
    return rcpp_result_gen;
END_RCPP
}
// guessBondsAtoms
Rcpp::DataFrame guessBondsAtoms(const Rcpp::S4& x, const Rcpp::Nullable< Rcpp::NumericVector > radius, const Rcpp::NumericVector& safety, const Rcpp::LogicalVector& usePBC);
RcppExport SEXP _Atoms_guessBondsAtoms(SEXP xSEXP, SEXP radiusSEXP, SEXP safetySEXP, SEXP usePBCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable< Rcpp::NumericVector > >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type safety(safetySEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type usePBC(usePBCSEXP);
    rcpp_result_gen = Rcpp::wrap(guessBondsAtoms(x, radius, safety, usePBC));
    return rcpp_result_gen;
END_RCPP
}
// guessDihedralsAtoms
Rcpp::DataFrame guessDihedralsAtoms(const Rcpp::S4& x);
RcppExport SEXP _Atoms_guessDihedralsAtoms(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(guessDihedralsAtoms(x));
    return rcpp_result_gen;
END_RCPP
}
// guessImpropersAtoms
Rcpp::DataFrame guessImpropersAtoms(const Rcpp::S4& x);
RcppExport SEXP _Atoms_guessImpropersAtoms(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(guessImpropersAtoms(x));
    return rcpp_result_gen;
END_RCPP
}
// joinAtoms
void joinAtoms(Rcpp::S4& x, const Rcpp::LogicalVector& multi);
RcppExport SEXP _Atoms_joinAtoms(SEXP xSEXP, SEXP multiSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type multi(multiSEXP);
    joinAtoms(x, multi);
    return R_NilValue;
END_RCPP
}
// measureBonds
Rcpp::NumericMatrix measureBonds(const Rcpp::S4& x, const Rcpp::IntegerVector& atm1, const Rcpp::IntegerVector& atm2);
RcppExport SEXP _Atoms_measureBonds(SEXP xSEXP, SEXP atm1SEXP, SEXP atm2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm1(atm1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm2(atm2SEXP);
    rcpp_result_gen = Rcpp::wrap(measureBonds(x, atm1, atm2));
    return rcpp_result_gen;
END_RCPP
}
// measureAngles
Rcpp::NumericMatrix measureAngles(const Rcpp::S4& x, const Rcpp::IntegerVector& atm1, const Rcpp::IntegerVector& atm2, const Rcpp::IntegerVector& atm3);
RcppExport SEXP _Atoms_measureAngles(SEXP xSEXP, SEXP atm1SEXP, SEXP atm2SEXP, SEXP atm3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm1(atm1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm2(atm2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm3(atm3SEXP);
    rcpp_result_gen = Rcpp::wrap(measureAngles(x, atm1, atm2, atm3));
    return rcpp_result_gen;
END_RCPP
}
// measureDihedrals
Rcpp::NumericMatrix measureDihedrals(const Rcpp::S4& x, const Rcpp::IntegerVector& atm1, const Rcpp::IntegerVector& atm2, const Rcpp::IntegerVector& atm3, const Rcpp::IntegerVector& atm4);
RcppExport SEXP _Atoms_measureDihedrals(SEXP xSEXP, SEXP atm1SEXP, SEXP atm2SEXP, SEXP atm3SEXP, SEXP atm4SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm1(atm1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm2(atm2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm3(atm3SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm4(atm4SEXP);
    rcpp_result_gen = Rcpp::wrap(measureDihedrals(x, atm1, atm2, atm3, atm4));
    return rcpp_result_gen;
END_RCPP
}
// planNormalsAtoms
Rcpp::S4 planNormalsAtoms(const Rcpp::S4& x, const Rcpp::IntegerVector& atm1, const Rcpp::IntegerVector& atm2, const Rcpp::IntegerVector& atm3);
RcppExport SEXP _Atoms_planNormalsAtoms(SEXP xSEXP, SEXP atm1SEXP, SEXP atm2SEXP, SEXP atm3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm1(atm1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm2(atm2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type atm3(atm3SEXP);
    rcpp_result_gen = Rcpp::wrap(planNormalsAtoms(x, atm1, atm2, atm3));
    return rcpp_result_gen;
END_RCPP
}
// rdfAtoms
Rcpp::DataFrame rdfAtoms(const Rcpp::S4& x, const Rcpp::IntegerVector& sel1, const Rcpp::IntegerVector& sel2, const Rcpp::IntegerVector& resnumb, const Rcpp::NumericVector& cutoff, const Rcpp::NumericVector& interval, const std::vector< bool >& type);
RcppExport SEXP _Atoms_rdfAtoms(SEXP xSEXP, SEXP sel1SEXP, SEXP sel2SEXP, SEXP resnumbSEXP, SEXP cutoffSEXP, SEXP intervalSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type sel1(sel1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type sel2(sel2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type resnumb(resnumbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< const std::vector< bool >& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(rdfAtoms(x, sel1, sel2, resnumb, cutoff, interval, type));
    return rcpp_result_gen;
END_RCPP
}
// rmsdAtoms
Rcpp::NumericVector rmsdAtoms(const Rcpp::S4& x, const Rcpp::S4& y, const Rcpp::LogicalVector& usePBC);
RcppExport SEXP _Atoms_rmsdAtoms(SEXP xSEXP, SEXP ySEXP, SEXP usePBCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type usePBC(usePBCSEXP);
    rcpp_result_gen = Rcpp::wrap(rmsdAtoms(x, y, usePBC));
    return rcpp_result_gen;
END_RCPP
}
// vectOrientation
Rcpp::NumericVector vectOrientation(const Rcpp::S4& x, const Rcpp::DataFrame& V1);
RcppExport SEXP _Atoms_vectOrientation(SEXP xSEXP, SEXP V1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V1(V1SEXP);
    rcpp_result_gen = Rcpp::wrap(vectOrientation(x, V1));
    return rcpp_result_gen;
END_RCPP
}
// vectOrientationAxis
Rcpp::NumericVector vectOrientationAxis(const Rcpp::S4& x, const Rcpp::DataFrame& V1, const Rcpp::NumericVector& V2);
RcppExport SEXP _Atoms_vectOrientationAxis(SEXP xSEXP, SEXP V1SEXP, SEXP V2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type V2(V2SEXP);
    rcpp_result_gen = Rcpp::wrap(vectOrientationAxis(x, V1, V2));
    return rcpp_result_gen;
END_RCPP
}
// vectOrientationPairwise
Rcpp::NumericVector vectOrientationPairwise(const Rcpp::S4& x, const Rcpp::DataFrame& V1, const Rcpp::DataFrame& V2);
RcppExport SEXP _Atoms_vectOrientationPairwise(SEXP xSEXP, SEXP V1SEXP, SEXP V2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V2(V2SEXP);
    rcpp_result_gen = Rcpp::wrap(vectOrientationPairwise(x, V1, V2));
    return rcpp_result_gen;
END_RCPP
}
// vectOrientationNotPairwise
Rcpp::NumericVector vectOrientationNotPairwise(const Rcpp::S4& x, const Rcpp::DataFrame& V1, const Rcpp::DataFrame& V2);
RcppExport SEXP _Atoms_vectOrientationNotPairwise(SEXP xSEXP, SEXP V1SEXP, SEXP V2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type V2(V2SEXP);
    rcpp_result_gen = Rcpp::wrap(vectOrientationNotPairwise(x, V1, V2));
    return rcpp_result_gen;
END_RCPP
}
// wrapAtoms
void wrapAtoms(Rcpp::S4& x, const Rcpp::IntegerVector& resnumb, const Rcpp::LogicalVector& multi);
RcppExport SEXP _Atoms_wrapAtoms(SEXP xSEXP, SEXP resnumbSEXP, SEXP multiSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type resnumb(resnumbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type multi(multiSEXP);
    wrapAtoms(x, resnumb, multi);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Atoms_LAMMPSParser", (DL_FUNC) &_Atoms_LAMMPSParser, 5},
    {"_Atoms_MOL2Parser", (DL_FUNC) &_Atoms_MOL2Parser, 5},
    {"_Atoms_PDBParser", (DL_FUNC) &_Atoms_PDBParser, 5},
    {"_Atoms_PDBWriter", (DL_FUNC) &_Atoms_PDBWriter, 2},
    {"_Atoms_RESParser", (DL_FUNC) &_Atoms_RESParser, 5},
    {"_Atoms_RESWriter", (DL_FUNC) &_Atoms_RESWriter, 2},
    {"_Atoms_XYZParser", (DL_FUNC) &_Atoms_XYZParser, 5},
    {"_Atoms_XYZWriter", (DL_FUNC) &_Atoms_XYZWriter, 2},
    {"_Atoms_comAtoms", (DL_FUNC) &_Atoms_comAtoms, 3},
    {"_Atoms_freeVolumeAtoms", (DL_FUNC) &_Atoms_freeVolumeAtoms, 2},
    {"_Atoms_freeVolume2Atoms", (DL_FUNC) &_Atoms_freeVolume2Atoms, 2},
    {"_Atoms_guessAnglesAtoms", (DL_FUNC) &_Atoms_guessAnglesAtoms, 1},
    {"_Atoms_guessBondsAtoms", (DL_FUNC) &_Atoms_guessBondsAtoms, 4},
    {"_Atoms_guessDihedralsAtoms", (DL_FUNC) &_Atoms_guessDihedralsAtoms, 1},
    {"_Atoms_guessImpropersAtoms", (DL_FUNC) &_Atoms_guessImpropersAtoms, 1},
    {"_Atoms_joinAtoms", (DL_FUNC) &_Atoms_joinAtoms, 2},
    {"_Atoms_measureBonds", (DL_FUNC) &_Atoms_measureBonds, 3},
    {"_Atoms_measureAngles", (DL_FUNC) &_Atoms_measureAngles, 4},
    {"_Atoms_measureDihedrals", (DL_FUNC) &_Atoms_measureDihedrals, 5},
    {"_Atoms_planNormalsAtoms", (DL_FUNC) &_Atoms_planNormalsAtoms, 4},
    {"_Atoms_rdfAtoms", (DL_FUNC) &_Atoms_rdfAtoms, 7},
    {"_Atoms_rmsdAtoms", (DL_FUNC) &_Atoms_rmsdAtoms, 3},
    {"_Atoms_vectOrientation", (DL_FUNC) &_Atoms_vectOrientation, 2},
    {"_Atoms_vectOrientationAxis", (DL_FUNC) &_Atoms_vectOrientationAxis, 3},
    {"_Atoms_vectOrientationPairwise", (DL_FUNC) &_Atoms_vectOrientationPairwise, 3},
    {"_Atoms_vectOrientationNotPairwise", (DL_FUNC) &_Atoms_vectOrientationNotPairwise, 3},
    {"_Atoms_wrapAtoms", (DL_FUNC) &_Atoms_wrapAtoms, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Atoms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}