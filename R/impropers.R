#' #' impropers
#' #' 
#' #' 
#' #' @name impropers
#' #' @export
#' impropers <- function(x, ...)
#'   UseMethod("impropers")
#' 
#' #' @rdname impropers
#' #' @export
#' impropers.Atoms <- function(x, atm1, atm2, atm3, atm4, atmprop = x$atmtype, pairwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2) && missing(atm3) && missing(atm4)){
#'     x@impropers
#'   } else if(!missing(atm1) && !missing(atm2) && !missing(atm3) && !missing(atm4)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(pairwise) || length(pairwise) != 1L)
#'       stop("'pairwise' must be a logical vector of length 1")
#'     if(pairwise){
#'       if(length(atm1) != length(atm2) ||
#'          length(atm1) != length(atm3) ||
#'          length(atm1) != length(atm4) )
#'         stop("'atm1', 'atm2', 'atm3' and 'atm4' ",
#'              "must have the same length when pairwise == TRUE")
#'       I <- paste(
#'         atmprop[x@impropers$atm1],
#'         atmprop[x@impropers$atm2],
#'         atmprop[x@impropers$atm3],
#'         atmprop[x@impropers$atm4])
#'       R1 <- paste(atm1, atm2, atm3, atm4)
#'       R2 <- paste(atm4, atm3, atm2, atm1)
#'       x@impropers[I %in% R1 | I %in% R2, ]
#'     } else {
#'       x@impropers[
#'         (atmprop[x@impropers$atm1] %in% atm1 &
#'          atmprop[x@impropers$atm2] %in% atm2 &
#'          atmprop[x@impropers$atm3] %in% atm3 &
#'          atmprop[x@impropers$atm4] %in% atm4) |
#'         (atmprop[x@impropers$atm4] %in% atm1 &
#'          atmprop[x@impropers$atm3] %in% atm2 &
#'          atmprop[x@impropers$atm2] %in% atm3 &
#'          atmprop[x@impropers$atm1] %in% atm4) , ]
#'     }
#'   } else {
#'     stop("'atm1', 'atm2', 'atm3' and 'atm4' must be specified")
#'   }
#' }
#' 
#' #' @rdname impropers
#' #' @export
#' "impropers<-" <- function(x, value)
#'   UseMethod("impropers<-")
#' 
#' #' @rdname impropers
#' #' @export
#' "impropers<-.Atoms" <- function(x, value){
#'   if(!is.data.frame(value))
#'     stop("'value' must be a data.frame")
#'   if(is.null(value$atm1))
#'     stop("'value' must have a 'atm1' column")
#'   if(is.null(value$atm2))
#'     stop("'value' must have a 'atm2' column")
#'   if(is.null(value$atm3))
#'     stop("'value' must have a 'atm3' column")
#'   if(is.null(value$atm4))
#'     stop("'value' must have a 'atm4' column")
#'   if(!is.integer(value$atm1))
#'     stop("'value$atm1' must be an integer vector")
#'   if(!is.integer(value$atm2))
#'     stop("'value$atm2' must be an integer vector")
#'   if(!is.integer(value$atm3))
#'     stop("'value$atm3' must be an integer vector")
#'   if(!is.integer(value$atm4))
#'     stop("'value$atm4' must be an integer vector")
#'   if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
#'      any(value$atm2 < 1) || any(value$atm2 > natom(x)) ||
#'      any(value$atm3 < 1) || any(value$atm3 > natom(x)) ||
#'      any(value$atm4 < 1) || any(value$atm4 > natom(x)) )
#'     stop("'value' contains atom indexes out of range 1:", natom(x))
#'   x@impropers <- data.frame(
#'     atm1 = value$atm1,
#'     atm2 = value$atm2,
#'     atm3 = value$atm3,
#'     atm4 = value$atm4)
#'   x@impropers[x@impropers$atm3 < x@impropers$atm2, ] <-
#'     x@impropers[x@impropers$atm3 < x@impropers$atm2, 4:1]
#'   x@impropers <- x@impropers[order(x@impropers$atm2, x@impropers$atm3), ]
#'   return(x)
#' }
