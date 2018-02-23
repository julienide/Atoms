#' planNormals
#' 
#' planNormals
#' 
#' @include checkSel.R
#' @name planNormals
#' @export
planNormals <- function(x, ...)
  UseMethod("planNormals")

#' @rdname planNormals
#' @export
planNormals.Atoms <- function(x, atm1, atm2 = NULL, atm3 = NULL, ...){
  if(!is.null(dim(atm1))){
    if(ncol(atm1) < 3L)
      stop("'atm1' must have at least three columns")
    atm3 <- atm1[, 3L]
    atm2 <- atm1[, 2L]
    atm1 <- atm1[, 1L]
  }
  natm <- natom(x)
  atm1 <- checkSel(atm1, natm)
  atm2 <- checkSel(atm2, natm)
  atm3 <- checkSel(atm3, natm)
  if(length(atm1) != length(atm2) ||
     length(atm1) != length(atm3))
    stop("'atm1', 'atm2' and 'atm3' must have the same length")
  ill.angles <-
    any(atm1 == atm2) ||
    any(atm2 == atm3) ||
    any(atm3 == atm1)
  if(ill.angles)
    warning("Ill defined plans")
  .planNormalsAtoms(x, atm1, atm2, atm3)
}
