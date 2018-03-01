#' #' Get or Set Bond Angles
#' #' 
#' #' Get or Set the indices of the atoms forming bond angles.
#' #' 
#' #' \code{angles} gets the indices of the atoms forming the bond angles specified
#' #' by \code{atm1}, \code{atm2} and \code{atm3}.
#' #' 
#' #' \code{angles<-} sets the indices of the atoms forming the bond angles.
#' #' 
#' #' 
#' #' @param x an R object
#' #' @param atm1,atm2,atm3 atomic properties of atoms forming bond angles.
#' #' @param atmprop a vector of atomic properties used for matching \code{atm1}, 
#' #'   \code{atm2}, \code{atm3}.
#' #' @param triplewise a logical value. Whether to consider \code{atm1}, 
#' #'   \code{atm2}, \code{atm3} as triplet or not.
#' #' @param value a three-column data.frame containing atom indices used to set 
#' #'   bond angles. Columns must be named \code{atm1}, \code{atm2}, \code{atm3}.
#' #'   
#' #' @return \code{angles} returns a three-column data.frame containing the
#' #'   indices of the atoms forming specific bond angles.
#' #'   
#' #'   \code{angles<-} returns an object of the same class as \code{x} with updated bond angles.
#' #' 
#' #' @name angles
#' #' @export
#' angles <- function(x, ...)
#'   UseMethod("angles")
#' 
#' #' @rdname angles
#' #' @export
#' angles.Atoms <- function(x, atm1, atm2, atm3, atmprop = x$atmtype, triplewise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2) && missing(atm3)){
#'     x@angles
#'   } else if(!missing(atm1) && !missing(atm2) && !missing(atm3)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(triplewise) || length(triplewise) != 1L)
#'       stop("'triplewise' must be a logical vector of length 1")
#'     if(triplewise){
#'       if(length(atm1) != length(atm2) ||
#'          length(atm1) != length(atm3) )
#'         stop("'atm1', 'atm2' and 'atm3' ",
#'              "must have the same length when triplewise == TRUE")
#'       A <- paste(
#'         atmprop[x@angles$atm1],
#'         atmprop[x@angles$atm2],
#'         atmprop[x@angles$atm3])
#'       R1 <- paste(atm1, atm2, atm3)
#'       R2 <- paste(atm3, atm2, atm1)
#'       x@angles[A %in% R1 | A %in% R2, ]
#'     } else {
#'       x@angles[
#'         (atmprop[x@angles$atm1] %in% atm1 &
#'          atmprop[x@angles$atm2] %in% atm2 &
#'          atmprop[x@angles$atm3] %in% atm3) |
#'         (atmprop[x@angles$atm3] %in% atm1 &
#'          atmprop[x@angles$atm2] %in% atm2 &
#'          atmprop[x@angles$atm1] %in% atm3) , ]
#'     }
#'   } else {
#'     stop("'atm1', 'atm2' and 'atm3' must be specified")
#'   }
#' }
#' 
#' #' @rdname angles
#' #' @export
#' "angles<-" <- function(x, value)
#'   UseMethod("angles<-")
#' 
#' #' @rdname angles
#' #' @export
#' "angles<-.Atoms" <- function(x, value){
#'   if(!is.data.frame(value))
#'     stop("'value' must be a data.frame")
#'   if(is.null(value$atm1))
#'     stop("'value' must have a 'atm1' column")
#'   if(is.null(value$atm2))
#'     stop("'value' must have a 'atm2' column")
#'   if(is.null(value$atm3))
#'     stop("'value' must have a 'atm3' column")
#'   if(!is.integer(value$atm1))
#'     stop("'value$atm1' must be an integer vector")
#'   if(!is.integer(value$atm2))
#'     stop("'value$atm2' must be an integer vector")
#'   if(!is.integer(value$atm3))
#'     stop("'value$atm3' must be an integer vector")
#'   if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
#'      any(value$atm2 < 1) || any(value$atm2 > natom(x)) ||
#'      any(value$atm3 < 1) || any(value$atm3 > natom(x)) )
#'     stop("'value' contains atom indexes out of range 1:", natom(x))
#'   x@angles <- data.frame(
#'     atm1 = value$atm1,
#'     atm2 = value$atm2,
#'     atm3 = value$atm3)
#'   x@angles[x@angles$atm3 < x@angles$atm1, ] <-
#'     x@angles[x@angles$atm3 < x@angles$atm1, 3:1]
#'   x@angles <- x@angles[order(x@angles$atm2), ]
#'   return(x)
#' }
#