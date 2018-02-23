#' Number Of Frames
#' 
#' Get the number of frames stored inside an objet.
#' 
#' The S4 method for class \sQuote{Atoms} count the number of columns of slot
#' \code{a}.
#' 
#' @param x an R object
#'   
#' @return an integer. The number of frames.
#' 
#' @seealso \code{\link{Atoms-class}}
#' 
#' @name nframe
#' @export
setGeneric(name = "nframe", def = function(x) standardGeneric("nframe"))

#' @export
setMethod(
  f = "nframe",
  signature = "Atoms",
  definition = function(x){
    ncol(x@a)
  }
)
