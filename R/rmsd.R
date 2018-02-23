#' Root Mean Square Deviation
#' 
#' Compute the root mean square deviation between a structure and a reference
#' structure.
#' 
#' \code{y} can be an object of class \sQuote{Atoms} or an integer indicating a
#' frame of \code{x} used as a reference.
#' 
#' \code{x} and \code{y} must have the same lattice vectors and PBC.
#' 
#' 
#' @param x an R object for which to compute the rmsd.
#' @param y the reference used to compute the rmsd. See details.
#' @param usePBC a logical value. Whether to use PBC to compute inter-atomic
#'   distances.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @name rmsd
#' @export
rmsd <- function(x, y, ...)
  UseMethod("rmsd")

#' @rdname rmsd
#' @export
rmsd.Atoms <- function(x, y, usePBC = TRUE){
  if(!is.logical(usePBC) || length(usePBC) != 1L)
    stop("'usePBC' must be a logical vector of length 1")
  if(is.numeric(y))
    y <- as.integer(y)
  if(is.integer(y)){
    if(length(y) != 1L)
      stop("'y' must be an integer vector of length 1")
    if(y < 1 || y > nframe(x))
      stop("'y' is out of range [1, ", nframe(x),"]")
    y <- x[, y]
  } else {
    if(!inherits(y, "Atoms"))
      stop("'y' must be an object of class 'Atoms'")
    if(nframe(y) != 1L)
      stop("'y' must contain a single frame")
    if(usePBC){
      if(!all(pbc(x) == pbc(y)))
        stop("'x' and 'y' must have the same PBC")
      if(!all(cell(x) == cell(y)))
        stop("'x' and 'y' must have the same lattice vectors to use PBC")
    }
  }
  if(natom(x) != natom(y))
    stop("'x' and 'y' must contain the same number of atom")
  
  .rmsdAtoms(x, y, usePBC)
}
