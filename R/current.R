#' Current Frame
#' 
#' Get the index of the current frame.
#' 
#' @param x an R object.
#' @param value an integer. Index of the current frame.
#'   
#' @return \code{current} returns an integer.
#'   
#'   \code{current<-} returns an object of the same class as \code{x} with the
#'   current frame index updated.
#'   
#' @export
setGeneric(name = "current", def = function(x) standardGeneric("current"))

#' @export
setMethod(
  f = "current",
  signature = "Atoms",
  definition = function(x){
    x@current
  }
)

#' @export
setGeneric(name = "current<-", def = function(x, value) standardGeneric("current<-"))

#' @export
setMethod(
  f = "current<-",
  signature = "Atoms",
  definition = function(x, value){
    if(!is.numeric(value) || length(value) != 1L)
      stop("'value' must be an integer vector of length 1")
    value <- as.integer(value)
    if(value < 1L || value > nframe(x))
      stop("'value' is out of frame range")
    x@current <- value
    return(x)
  }
)
