#' #' Extract or Replace Topology Components
#' #'
#' #' Extract or Replace bonds, bond angles, dihedral angles and improper angles.
#' #'
#' #' @param x an R object
#' #' @param atm1,atm2,atm3,atm4 atomic properties of atoms forming bonds, bond
#' #'   angles, dihedral angles or improper angles.
#' #' @param atmprop a vector of atomic properties used for matching \code{atm1},
#' #'   \code{atm2}, \code{atm3} and \code{atm4}.
#' #' @param nwise a logical value. Whether to consider \code{atm1} and
#' #'   \code{atm2}/\code{atm3}/\code{atm4} as pair/triplet/quadruplet or not.
#' #' @param value a data.frame containing by column \code{atm1}, \code{atm2},
#' #'   \code{atm3} and \code{atm4} indices used to set bonds, bond angles,
#' #'   dihedral angles or improper angles. Columns must be named \code{atm1},
#' #'   \code{atm2}, \code{atm3} and \code{atm4}.
#' #'
#' #' @return \code{bonds}, \code{angles}, \code{dihedrals} and \code{impropers}
#' #'   return a data.frame containing the indices of the atoms forming specific
#' #'   bonds, bond angles, dihedral angles, improper angles.
#' #'
#' #'   \code{bonds<-}, \code{angles<-}, \code{dihedrals<-} and \code{impropers<-}
#' #'   return an object of the same class as \code{x} with updated topology.
#' #'
#' #'
#' #'
#' #' @name topologyGettersAndSetters
#' #' @export
#' bonds <- function(x, ...)
#'   UseMethod("bonds")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' angles <- function(x, ...)
#'   UseMethod("angles")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' dihedrals <- function(x, ...)
#'   UseMethod("dihedrals")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' impropers <- function(x, ...)
#'   UseMethod("impropers")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' "bonds<-" <- function(x, value)
#'   UseMethod("bonds<-")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' "angles<-" <- function(x, value)
#'   UseMethod("angles<-")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' "dihedrals<-" <- function(x, value)
#'   UseMethod("dihedrals<-")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' "impropers<-" <- function(x, value)
#'   UseMethod("impropers<-")
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' bonds.Atoms <- function(x, atm1, atm2, atmprop = x$atmtype, nwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2)){
#'     x@bonds
#'   } else if(!missing(atm1) && !missing(atm2)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(nwise) || length(nwise) != 1L)
#'       stop("'nwise' must be a logical vector of length 1")
#'     if(nwise){
#'       if(length(atm1) != length(atm2))
#'         stop("'atm1' and 'atm2' ",
#'              "must have the same length when nwise == TRUE")
#'       B <- paste(
#'         atmprop[x@bonds$atm1],
#'         atmprop[x@bonds$atm2])
#'       R1 <- paste(atm1, atm2)
#'       R2 <- paste(atm2, atm1)
#'       B <- x@bonds[B %in% R1 | B %in% R2, ]
#'     } else {
#'       B <- x@bonds[
#'         (atmprop[x@bonds$atm1] %in% atm1 &
#'            atmprop[x@bonds$atm2] %in% atm2) |
#'           (atmprop[x@bonds$atm2] %in% atm1 &
#'              atmprop[x@bonds$atm1] %in% atm2), ]
#'     }
#'     B[!(atmprop[B$atm1] %in% atm1), ] <- B[!(atmprop[B$atm1] %in% atm1), 2:1]
#'     return(B)
#'   } else {
#'     stop("'atm1' and 'atm2' must be specified")
#'   }
#' }
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' angles.Atoms <- function(x, atm1, atm2, atm3, atmprop = x$atmtype, nwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2) && missing(atm3)){
#'     x@angles
#'   } else if(!missing(atm1) && !missing(atm2) && !missing(atm3)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(nwise) || length(nwise) != 1L)
#'       stop("'nwise' must be a logical vector of length 1")
#'     if(nwise){
#'       if(length(atm1) != length(atm2) ||
#'          length(atm1) != length(atm3) )
#'         stop("'atm1', 'atm2' and 'atm3' ",
#'              "must have the same length when nwise == TRUE")
#'       A <- paste(
#'         atmprop[x@angles$atm1],
#'         atmprop[x@angles$atm2],
#'         atmprop[x@angles$atm3])
#'       R1 <- paste(atm1, atm2, atm3)
#'       R2 <- paste(atm3, atm2, atm1)
#'       A <- x@angles[A %in% R1 | A %in% R2, ]
#'     } else {
#'       A <- x@angles[
#'         (atmprop[x@angles$atm1] %in% atm1 &
#'            atmprop[x@angles$atm2] %in% atm2 &
#'            atmprop[x@angles$atm3] %in% atm3) |
#'           (atmprop[x@angles$atm3] %in% atm1 &
#'              atmprop[x@angles$atm2] %in% atm2 &
#'              atmprop[x@angles$atm1] %in% atm3) , ]
#'     }
#'     A[!(atmprop[A$atm1] %in% atm1), ] <- A[!(atmprop[A$atm1] %in% atm1), 3:1]
#'     return(A)
#'   } else {
#'     stop("'atm1', 'atm2' and 'atm3' must be specified")
#'   }
#' }
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' dihedrals.Atoms <- function(x, atm1, atm2, atm3, atm4, atmprop = x$atmtype, nwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2) && missing(atm3) && missing(atm4)){
#'     x@dihedrals
#'   } else if(!missing(atm1) && !missing(atm2) && !missing(atm3) && !missing(atm4)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(nwise) || length(nwise) != 1L)
#'       stop("'nwise' must be a logical vector of length 1")
#'     if(nwise){
#'       if(length(atm1) != length(atm2) ||
#'          length(atm1) != length(atm3) ||
#'          length(atm1) != length(atm4) )
#'         stop("'atm1', 'atm2', 'atm3' and 'atm4' ",
#'              "must have the same length when nwise == TRUE")
#'       D <- paste(
#'         atmprop[x@dihedrals$atm1],
#'         atmprop[x@dihedrals$atm2],
#'         atmprop[x@dihedrals$atm3],
#'         atmprop[x@dihedrals$atm4])
#'       R1 <- paste(atm1, atm2, atm3, atm4)
#'       R2 <- paste(atm4, atm3, atm2, atm1)
#'       D <- x@dihedrals[D %in% R1 | D %in% R2, ]
#'     } else {
#'       D <- x@dihedrals[
#'         (atmprop[x@dihedrals$atm1] %in% atm1 &
#'            atmprop[x@dihedrals$atm2] %in% atm2 &
#'            atmprop[x@dihedrals$atm3] %in% atm3 &
#'            atmprop[x@dihedrals$atm4] %in% atm4) |
#'           (atmprop[x@dihedrals$atm4] %in% atm1 &
#'              atmprop[x@dihedrals$atm3] %in% atm2 &
#'              atmprop[x@dihedrals$atm2] %in% atm3 &
#'              atmprop[x@dihedrals$atm1] %in% atm4) , ]
#'     }
#'     D[!(atmprop[D$atm1] %in% atm1), ] <- D[!(atmprop[D$atm1] %in% atm1), 4:1]
#'     return(D)
#'   } else {
#'     stop("'atm1', 'atm2', 'atm3' and 'atm4' must be specified")
#'   }
#' }
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' impropers.Atoms <- function(x, atm1, atm2, atm3, atm4, atmprop = x$atmtype, nwise = TRUE, ...){
#'   if(missing(atm1) && missing(atm2) && missing(atm3) && missing(atm4)){
#'     x@impropers
#'   } else if(!missing(atm1) && !missing(atm2) && !missing(atm3) && !missing(atm4)){
#'     if(!is.vector(atmprop))
#'       stop("'atmprop' must be a vector")
#'     if(length(atmprop) != natom(x))
#'       stop("'atmprop' must be of length ", natom(x))
#'     if(!is.logical(nwise) || length(nwise) != 1L)
#'       stop("'nwise' must be a logical vector of length 1")
#'     if(nwise){
#'       if(length(atm1) != length(atm2) ||
#'          length(atm1) != length(atm3) ||
#'          length(atm1) != length(atm4) )
#'         stop("'atm1', 'atm2', 'atm3' and 'atm4' ",
#'              "must have the same length when nwise == TRUE")
#'       I <- paste(
#'         atmprop[x@impropers$atm1],
#'         atmprop[x@impropers$atm2],
#'         atmprop[x@impropers$atm3],
#'         atmprop[x@impropers$atm4])
#'       R1 <- paste(atm1, atm2, atm3, atm4)
#'       R2 <- paste(atm4, atm3, atm2, atm1)
#'       I <- x@impropers[I %in% R1 | I %in% R2, ]
#'     } else {
#'       I <- x@impropers[
#'         (atmprop[x@impropers$atm1] %in% atm1 &
#'            atmprop[x@impropers$atm2] %in% atm2 &
#'            atmprop[x@impropers$atm3] %in% atm3 &
#'            atmprop[x@impropers$atm4] %in% atm4) |
#'           (atmprop[x@impropers$atm4] %in% atm1 &
#'              atmprop[x@impropers$atm3] %in% atm2 &
#'              atmprop[x@impropers$atm2] %in% atm3 &
#'              atmprop[x@impropers$atm1] %in% atm4) , ]
#'     }
#'     I[!(atmprop[I$atm1] %in% atm1), ] <- I[!(atmprop[I$atm1] %in% atm1), 4:1]
#'     return(I)
#'   } else {
#'     stop("'atm1', 'atm2', 'atm3' and 'atm4' must be specified")
#'   }
#' }
#' 
#' #' @rdname topologyGettersAndSetters
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
#'   rownames(x@bonds) <- NULL
#'   return(x)
#' }
#' 
#' #' @rdname topologyGettersAndSetters
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
#'   rownames(x@angles) <- NULL
#'   return(x)
#' }
#' 
#' #' @rdname topologyGettersAndSetters
#' #' @export
#' "dihedrals<-.Atoms" <- function(x, value){
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
#'   x@dihedrals <- data.frame(
#'     atm1 = value$atm1,
#'     atm2 = value$atm2,
#'     atm3 = value$atm3,
#'     atm4 = value$atm4)
#'   x@dihedrals[x@dihedrals$atm3 < x@dihedrals$atm2, ] <-
#'     x@dihedrals[x@dihedrals$atm3 < x@dihedrals$atm2, 4:1]
#'   x@dihedrals <- x@dihedrals[order(x@dihedrals$atm2, x@dihedrals$atm3), ]
#'   rownames(x@dihedrals) <- NULL
#'   return(x)
#' }
#' 
#' #' @rdname topologyGettersAndSetters
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
#'   rownames(x@impropers) <- NULL
#'   return(x)
#' }
