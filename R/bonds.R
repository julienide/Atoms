#' #' bonds
#' #' 
#' #' 
#' #' @name bonds
#' #' @export
#' bonds <- function(x, ...)
#'   UseMethod("bonds")
#' 
#' #' @rdname bonds
#' #' @export
#' bonds.Atoms <- function(x, atm1, atm2, atmprop = x$atmtype, pairwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2)){
#'     x@bonds
#'   } else if(!missing(atm1) && !missing(atm2)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     # if(length(atm1) != 1L)
#'     #   stop("'atm1' must be a vector of length 1")
#'     # if(length(atm2) != 1L)
#'     #   stop("'atm2' must be a vector of length 1")
#'     if(!is.logical(pairwise) || length(pairwise) != 1L)
#'       stop("'pairwise' must be a logical vector of length 1")
#'     if(pairwise){
#'       if(length(atm1) != length(atm2))
#'         stop("'atm1' and 'atm2' ",
#'              "must have the same length when pairwise == TRUE")
#'       B <- paste(
#'         atmprop[x@bonds$atm1],
#'         atmprop[x@bonds$atm2])
#'       R1 <- paste(atm1, atm2)
#'       R2 <- paste(atm2, atm1)
#'       x@bonds[B %in% R1 | B %in% R2, ]
#'     } else {
#'       x@bonds[
#'         (atmprop[x@bonds$atm1] %in% atm1 &
#'          atmprop[x@bonds$atm2] %in% atm2) |
#'         (atmprop[x@bonds$atm2] %in% atm1 &
#'          atmprop[x@bonds$atm1] %in% atm2), ]
#'     }
#'   } else {
#'     stop("'atm1' and 'atm2' must be specified")
#'   }
#' }
#' 
#' #' @rdname bonds
#' #' @export
#' "bonds<-" <- function(x, value)
#'   UseMethod("bonds<-")
#' 
#' #' @rdname bonds
#' #' @export
#' "bonds<-.Atoms" <- function(x, value){
#'   if(!is.data.frame(value))
#'     stop("'value' must be a data.frame")
#'   if(is.null(value$atm1))
#'     stop("'value' must have a 'atm1' column")
#'   if(is.null(value$atm2))
#'     stop("'value' must have a 'atm2' column")
#'   if(!is.integer(value$atm1))
#'     stop("'value$atm1' must be an integer vector")
#'   if(!is.integer(value$atm2))
#'     stop("'value$atm2' must be an integer vector")
#'   if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
#'      any(value$atm2 < 1) || any(value$atm2 > natom(x)) )
#'     stop("'value' contains atom indexes out of range 1:", natom(x))
#'   x@bonds <- data.frame(
#'     atm1 = value$atm1,
#'     atm2 = value$atm2)
#'   x@bonds[x@bonds$atm2 < x@bonds$atm1, ] <-
#'     x@bonds[x@bonds$atm2 < x@bonds$atm1, 2:1]
#'   x@bonds <- x@bonds[order(x@bonds$atm1, x@bonds$atm2), ]
#'   return(x)
#' }
