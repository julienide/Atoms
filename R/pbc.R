#' pbc
#'
#' @name pbc
#' @export
setGeneric(name = "pbc", def = function(x) standardGeneric("pbc"))

#' @rdname pbc
#' @export
setMethod(
  f = "pbc",
  signature = "Atoms",
  definition = function(x){
    x@pbc
  }
)

#' @rdname pbc
#' @export
setGeneric(name = "pbc<-", def = function(x, value) standardGeneric("pbc<-"))

#' @rdname pbc
#' @export
setMethod(
  f = "pbc<-",
  signature = "Atoms",
  definition = function(x, value){
    if(!is.logical(value) || !(length(value) %in% c(1L, 3L)))
      stop("'value' must be a logical vector of length 3")
    if(length(value) == 1L)
      value <- rep(value, 3L)
    x@pbc <- value
    # if(!x@pbc[1L])
    #   x@a <- c(1.0, 0.0, 0.0)
    # if(!x@pbc[2L])
    #   x@b <- c(0.0, 1.0, 0.0)
    # if(!x@pbc[3L])
    #   x@c <- c(0.0, 0.0, 1.0)
    return(x)
  }
)
