#' Atomic Recognition
#' 
#' Determine atomic symbols from an \sQuote{Atoms} object.
#' 
#' This function apply the \code{symb} method to the \code{x$atmtype} component 
#' of an \sQuote{Atoms} object.\cr If \code{x$atmtype} is \code{NULL}, 
#' \code{x$atmname} is used instead.\cr If \code{x$atmtype} contains empty 
#' character strings, they are replaced by their associated \code{x$atmname} 
#' values.
#' 
#' @return a character vector. Atomic symbols.
#'   
#' @param x an object of class \sQuote{Atoms}
#' @param na.as.dyummy a logical value. Whether to consider NA values as dummy
#'   atoms or not.
#'
#' @keywords manip
#' @name atomicRecognition
#' @export
symb.Atoms <- function(x, na.as.dummy = FALSE){
  if(!is.null(x$atmtype)){
    s <- x$atmtype
    if(!is.null(x$atmname)){
      s[!nzchar(s)] <- x$atmname[!nzchar(s)]
    }
  } else if(!si.null(x$atmname)) {
    s <- x$atmname
  } else {
    s <- rep("Xx", natom(x))
  }
  symb(s)
}
