#' Measure Bonds Angles or Dihedrals
#' 
#' 
#' @name measure
#' @include checkSel.R
#' @export
measure <- function(x, atm1, atm2, atm3, atm4)
  UseMethod("measure")

#' @rdname measure
#' @export
measure.Atoms <- function(x, atm1, atm2 = NULL, atm3 = NULL, atm4 = NULL){
  if(!is.null(dim(atm1))){
    if(ncol(atm1) < 2L)
      stop("'atm1' must have at least two columns")
    if(ncol(atm1) >= 4L)
      atm4 <- atm1[, 4L]
    if(ncol(atm1) >= 3L)
      atm3 <- atm1[, 3L]
    atm2 <- atm1[, 2L]
    atm1 <- atm1[, 1L]
  }
    
  natm <- natom(x)
  atm1 <- checkSel(atm1, natm)
  atm2 <- checkSel(atm2, natm)
  if(length(atm1) != length(atm2))
    stop("'atm1' and 'atm2' must have the same length")
  if(is.null(atm4)){
    if(is.null(atm3)){
      .measureBonds(x, atm1, atm2)
    } else {
      atm3 <- checkSel(atm3, natm)
      if(length(atm1) != length(atm3))
        stop("'atm1', 'atm2' and 'atm3' must have the same length")
      .measureAngles(x, atm1, atm2, atm3)
    }
  } else {
    atm3 <- checkSel(atm3, natm)
    atm4 <- checkSel(atm4, natm)
    if(length(atm1) != length(atm3) ||
       length(atm1) != length(atm4))
      stop("'atm1', 'atm2', 'atm3' and 'atm4' must have the same length")
    .measureDihedrals(x, atm1, atm2, atm3, atm4)
  }
}
