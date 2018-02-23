#' Number Of Atoms
#' 
#' Get the number of atoms stored inside an objet.
#' 
#' The S4 method for class \sQuote{Atoms} count the number of rows of slot
#' \code{atoms}.
#' 
#' @param x an R object
#'   
#' @return an integer. The number of atoms.
#' 
#' @seealso \code{\link{Atoms-class}}
#' 
#' @name natom
#' @export
setGeneric(name = "natom", def = function(x) standardGeneric("natom"))

#' @export
setMethod(
  f = "natom",
  signature = "Atoms",
  definition = function(x){
    nrow(x@atoms)
  }
)
