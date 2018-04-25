#' Extract or Replace Topology Components
#'
#' Extract or Replace bonds, bond angles, dihedral angles and improper angles.
#'
#' @param x an R object
#' @param atm1,atm2,atm3,atm4 atomic properties of atoms forming bonds, bond
#'   angles, dihedral angles or improper angles.
#' @param atmprop a vector of atomic properties used for matching \code{atm1},
#'   \code{atm2}, \code{atm3} and \code{atm4}.
#' @param nwise a logical value. Whether to consider \code{atm1} and
#'   \code{atm2}/\code{atm3}/\code{atm4} as pair/triplet/quadruplet or not.
#' @param value a data.frame containing by column \code{atm1}, \code{atm2},
#'   \code{atm3} and \code{atm4} indices used to set bonds, bond angles,
#'   dihedral angles or improper angles. Columns must be named \code{atm1},
#'   \code{atm2}, \code{atm3} and \code{atm4}.
#'
#' @return \code{bonds}, \code{angles}, \code{dihedrals} and \code{impropers}
#'   return a data.frame containing the indices of the atoms forming specific
#'   bonds, bond angles, dihedral angles, improper angles.
#'
#'   \code{bonds<-}, \code{angles<-}, \code{dihedrals<-} and \code{impropers<-}
#'   return an object of the same class as \code{x} with updated topology.
#'
#'
#'
#' @name topologyGettersAndSetters
#' @export
bonds <- function(x, ...)
  UseMethod("bonds")

#' @rdname topologyGettersAndSetters
#' @export
angles <- function(x, ...)
  UseMethod("angles")

#' @rdname topologyGettersAndSetters
#' @export
dihedrals <- function(x, ...)
  UseMethod("dihedrals")

#' @rdname topologyGettersAndSetters
#' @export
impropers <- function(x, ...)
  UseMethod("impropers")

#' @rdname topologyGettersAndSetters
#' @export
"bonds<-" <- function(x, value)
  UseMethod("bonds<-")

#' @rdname topologyGettersAndSetters
#' @export
"angles<-" <- function(x, value)
  UseMethod("angles<-")

#' @rdname topologyGettersAndSetters
#' @export
"dihedrals<-" <- function(x, value)
  UseMethod("dihedrals<-")

#' @rdname topologyGettersAndSetters
#' @export
"impropers<-" <- function(x, value)
  UseMethod("impropers<-")

#' @rdname topologyGettersAndSetters
#' @export
bonds.Atoms <- function(x, atm1 = NA, atm2 = NA,
                         atmprop = x$atmtype, nwise = TRUE, ...){
  if(!is.vector(atmprop))
    stop("'atmprop' must be a vector")
  if(length(atmprop) != natom(x))
    stop("'atmprop' must be of length ", natom(x))
  if(!is.logical(nwise) || length(nwise) != 1L)
    stop("'nwise' must be a logical vector of length 1")
  x <- x@bonds
  if(nwise){
    if(length(atm1) != length(atm2))
      stop("'atm1' and 'atm2' ",
           "must have the same length when nwise == TRUE")
    atm1[is.na(atm1)] <- ".*"
    atm2[is.na(atm2)] <- ".*"
    T1 <- paste("^", atm1, " ", atm2, "$", sep = "", collapse = "|")
    T2 <- paste("^", atm2, " ", atm1, "$", sep = "", collapse = "|")
    T1 <- grepl(T1, paste(atmprop[x$atm1], atmprop[x$atm2]))
    T2 <- grepl(T2, paste(atmprop[x$atm1], atmprop[x$atm2]))
  } else {
    if(anyNA(atm1) || anyNA(atm2))
      stop("'atm1' and 'atm2' can not contain NA values when nwise == FALSE")
    T1 <- atmprop[x$atm1] %in% atm1 & atmprop[x$atm2] %in% atm2
    T2 <- atmprop[x$atm2] %in% atm1 & atmprop[x$atm1] %in% atm2
  }
  x[!T1 & T2, ] <- x[!T1 & T2, 2:1]
  x[T1|T2, ]
}

#' @rdname topologyGettersAndSetters
#' @export
angles.Atoms <- function(x, atm1 = NA, atm2 = NA, atm3 = NA,
                         atmprop = x$atmtype, nwise = TRUE, ...){
  if(!is.vector(atmprop))
    stop("'atmprop' must be a vector")
  if(length(atmprop) != natom(x))
    stop("'atmprop' must be of length ", natom(x))
  if(!is.logical(nwise) || length(nwise) != 1L)
    stop("'nwise' must be a logical vector of length 1")
  x <- x@angles
  if(nwise){
    if(length(atm1) != length(atm2) ||
       length(atm1) != length(atm3) )
      stop("'atm1', 'atm2' and 'atm3' ",
           "must have the same length when nwise == TRUE")
    atm1[is.na(atm1)] <- ".*"
    atm2[is.na(atm2)] <- ".*"
    atm3[is.na(atm3)] <- ".*"
    T1 <- paste("^", atm1, " ", atm2, " ", atm3, "$", sep = "", collapse = "|")
    T2 <- paste("^", atm3, " ", atm2, " ", atm1, "$", sep = "", collapse = "|")
    T1 <- grepl(T1, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3]))
    T2 <- grepl(T2, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3]))
  } else {
    if(anyNA(atm1) || anyNA(atm2) || anyNA(atm3))
      stop("'atm1', 'atm2' and 'atm3' can not contain NA values when nwise == FALSE")
    T1 <- atmprop[x$atm1] %in% atm1 & atmprop[x$atm2] %in% atm2 & atmprop[x$atm3] %in% atm3
    T2 <- atmprop[x$atm3] %in% atm1 & atmprop[x$atm2] %in% atm2 & atmprop[x$atm1] %in% atm3
  }
  x[!T1 & T2, ] <- x[!T1 & T2, 3:1]
  x[T1|T2, ]
}

#' @rdname topologyGettersAndSetters
#' @export
dihedrals.Atoms <- function(x, atm1 = NA, atm2 = NA, atm3 = NA, atm4 = NA,
                            atmprop = x$atmtype, nwise = TRUE, ...){
  if(!is.vector(atmprop))
    stop("'atmprop' must be a vector")
  if(length(atmprop) != natom(x))
    stop("'atmprop' must be of length ", natom(x))
  if(!is.logical(nwise) || length(nwise) != 1L)
    stop("'nwise' must be a logical vector of length 1")
  x <- x@dihedrals
  if(nwise){
    if(length(atm1) != length(atm2) ||
       length(atm1) != length(atm3) ||
       length(atm1) != length(atm4) )
      stop("'atm1', 'atm2', 'atm3' and 'atm4' ",
           "must have the same length when nwise == TRUE")
    atm1[is.na(atm1)] <- ".*"
    atm2[is.na(atm2)] <- ".*"
    atm3[is.na(atm3)] <- ".*"
    atm4[is.na(atm4)] <- ".*"
    T1 <- paste("^", atm1, " ", atm2, " ", atm3, " ", atm4, "$", sep = "", collapse = "|")
    T2 <- paste("^", atm4, " ", atm3, " ", atm2, " ", atm1, "$", sep = "", collapse = "|")
    T1 <- grepl(T1, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3], atmprop[x$atm4]))
    T2 <- grepl(T2, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3], atmprop[x$atm4]))
  } else {
    if(anyNA(atm1) || anyNA(atm2) || anyNA(atm3) || anyNA(atm4))
      stop("'atm1', 'atm2', 'atm3' and 'atm4' can not contain NA values when nwise == FALSE")
    T1 <-
      atmprop[x$atm1] %in% atm1 &
      atmprop[x$atm2] %in% atm2 &
      atmprop[x$atm3] %in% atm3 &
      atmprop[x$atm4] %in% atm4
    T2 <-
      atmprop[x$atm4] %in% atm1 &
      atmprop[x$atm3] %in% atm2 &
      atmprop[x$atm2] %in% atm3 &
      atmprop[x$atm1] %in% atm4
  }
  x[!T1 & T2, ] <- x[!T1 & T2, 4:1]
  x[T1|T2, ]
}

#' @rdname topologyGettersAndSetters
#' @export
impropers.Atoms <- function(x, atm1 = NA, atm2 = NA, atm3 = NA, atm4 = NA,
                            atmprop = x$atmtype, nwise = TRUE, ...){
  if(!is.vector(atmprop))
    stop("'atmprop' must be a vector")
  if(length(atmprop) != natom(x))
    stop("'atmprop' must be of length ", natom(x))
  if(!is.logical(nwise) || length(nwise) != 1L)
    stop("'nwise' must be a logical vector of length 1")
  x <- x@impropers
  if(nwise){
    if(length(atm1) != length(atm2) ||
       length(atm1) != length(atm3) ||
       length(atm1) != length(atm4) )
      stop("'atm1', 'atm2', 'atm3' and 'atm4' ",
           "must have the same length when nwise == TRUE")
    atm1[is.na(atm1)] <- ".*"
    atm2[is.na(atm2)] <- ".*"
    atm3[is.na(atm3)] <- ".*"
    atm4[is.na(atm4)] <- ".*"
    T1 <- paste("^", atm1, " ", atm2, " ", atm3, " ", atm4, "$", sep = "", collapse = "|")
    T2 <- paste("^", atm4, " ", atm3, " ", atm2, " ", atm1, "$", sep = "", collapse = "|")
    T1 <- grepl(T1, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3], atmprop[x$atm4]))
    T2 <- grepl(T2, paste(atmprop[x$atm1], atmprop[x$atm2], atmprop[x$atm3], atmprop[x$atm4]))
  } else {
    if(anyNA(atm1) || anyNA(atm2) || anyNA(atm3) || anyNA(atm4))
      stop("'atm1', 'atm2', 'atm3' and 'atm4' can not contain NA values when nwise == FALSE")
    T1 <-
      atmprop[x$atm1] %in% atm1 &
      atmprop[x$atm2] %in% atm2 &
      atmprop[x$atm3] %in% atm3 &
      atmprop[x$atm4] %in% atm4
    T2 <-
      atmprop[x$atm4] %in% atm1 &
      atmprop[x$atm3] %in% atm2 &
      atmprop[x$atm2] %in% atm3 &
      atmprop[x$atm1] %in% atm4
  }
  x[!T1 & T2, ] <- x[!T1 & T2, 4:1]
  x[T1|T2, ]
}

#' @rdname topologyGettersAndSetters
#' @export
"bonds<-.Atoms" <- function(x, value){
  if(!is.data.frame(value))
    stop("'value' must be a data.frame")
  if(is.null(value$atm1))
    stop("'value' must have a 'atm1' column")
  if(is.null(value$atm2))
    stop("'value' must have a 'atm2' column")
  if(!is.integer(value$atm1))
    stop("'value$atm1' must be an integer vector")
  if(!is.integer(value$atm2))
    stop("'value$atm2' must be an integer vector")
  if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
     any(value$atm2 < 1) || any(value$atm2 > natom(x)) )
    stop("'value' contains atom indexes out of range 1:", natom(x))
  x@bonds <- data.frame(
    atm1 = value$atm1,
    atm2 = value$atm2)
  x@bonds[x@bonds$atm2 < x@bonds$atm1, ] <-
    x@bonds[x@bonds$atm2 < x@bonds$atm1, 2:1]
  x@bonds <- x@bonds[order(x@bonds$atm1, x@bonds$atm2), ]
  rownames(x@bonds) <- NULL
  return(x)
}

#' @rdname topologyGettersAndSetters
#' @export
"angles<-.Atoms" <- function(x, value){
  if(!is.data.frame(value))
    stop("'value' must be a data.frame")
  if(is.null(value$atm1))
    stop("'value' must have a 'atm1' column")
  if(is.null(value$atm2))
    stop("'value' must have a 'atm2' column")
  if(is.null(value$atm3))
    stop("'value' must have a 'atm3' column")
  if(!is.integer(value$atm1))
    stop("'value$atm1' must be an integer vector")
  if(!is.integer(value$atm2))
    stop("'value$atm2' must be an integer vector")
  if(!is.integer(value$atm3))
    stop("'value$atm3' must be an integer vector")
  if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
     any(value$atm2 < 1) || any(value$atm2 > natom(x)) ||
     any(value$atm3 < 1) || any(value$atm3 > natom(x)) )
    stop("'value' contains atom indexes out of range 1:", natom(x))
  x@angles <- data.frame(
    atm1 = value$atm1,
    atm2 = value$atm2,
    atm3 = value$atm3)
  x@angles[x@angles$atm3 < x@angles$atm1, ] <-
    x@angles[x@angles$atm3 < x@angles$atm1, 3:1]
  x@angles <- x@angles[order(x@angles$atm2), ]
  rownames(x@angles) <- NULL
  return(x)
}

#' @rdname topologyGettersAndSetters
#' @export
"dihedrals<-.Atoms" <- function(x, value){
  if(!is.data.frame(value))
    stop("'value' must be a data.frame")
  if(is.null(value$atm1))
    stop("'value' must have a 'atm1' column")
  if(is.null(value$atm2))
    stop("'value' must have a 'atm2' column")
  if(is.null(value$atm3))
    stop("'value' must have a 'atm3' column")
  if(is.null(value$atm4))
    stop("'value' must have a 'atm4' column")
  if(!is.integer(value$atm1))
    stop("'value$atm1' must be an integer vector")
  if(!is.integer(value$atm2))
    stop("'value$atm2' must be an integer vector")
  if(!is.integer(value$atm3))
    stop("'value$atm3' must be an integer vector")
  if(!is.integer(value$atm4))
    stop("'value$atm4' must be an integer vector")
  if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
     any(value$atm2 < 1) || any(value$atm2 > natom(x)) ||
     any(value$atm3 < 1) || any(value$atm3 > natom(x)) ||
     any(value$atm4 < 1) || any(value$atm4 > natom(x)) )
    stop("'value' contains atom indexes out of range 1:", natom(x))
  x@dihedrals <- data.frame(
    atm1 = value$atm1,
    atm2 = value$atm2,
    atm3 = value$atm3,
    atm4 = value$atm4)
  x@dihedrals[x@dihedrals$atm3 < x@dihedrals$atm2, ] <-
    x@dihedrals[x@dihedrals$atm3 < x@dihedrals$atm2, 4:1]
  x@dihedrals <- x@dihedrals[order(x@dihedrals$atm2, x@dihedrals$atm3), ]
  rownames(x@dihedrals) <- NULL
  return(x)
}

#' @rdname topologyGettersAndSetters
#' @export
"impropers<-.Atoms" <- function(x, value){
  if(!is.data.frame(value))
    stop("'value' must be a data.frame")
  if(is.null(value$atm1))
    stop("'value' must have a 'atm1' column")
  if(is.null(value$atm2))
    stop("'value' must have a 'atm2' column")
  if(is.null(value$atm3))
    stop("'value' must have a 'atm3' column")
  if(is.null(value$atm4))
    stop("'value' must have a 'atm4' column")
  if(!is.integer(value$atm1))
    stop("'value$atm1' must be an integer vector")
  if(!is.integer(value$atm2))
    stop("'value$atm2' must be an integer vector")
  if(!is.integer(value$atm3))
    stop("'value$atm3' must be an integer vector")
  if(!is.integer(value$atm4))
    stop("'value$atm4' must be an integer vector")
  if(any(value$atm1 < 1) || any(value$atm1 > natom(x)) ||
     any(value$atm2 < 1) || any(value$atm2 > natom(x)) ||
     any(value$atm3 < 1) || any(value$atm3 > natom(x)) ||
     any(value$atm4 < 1) || any(value$atm4 > natom(x)) )
    stop("'value' contains atom indexes out of range 1:", natom(x))
  x@impropers <- data.frame(
    atm1 = value$atm1,
    atm2 = value$atm2,
    atm3 = value$atm3,
    atm4 = value$atm4)
  # x@impropers[x@impropers$atm3 < x@impropers$atm2, ] <-
  #   x@impropers[x@impropers$atm3 < x@impropers$atm2, 4:1]
  x@impropers <- x@impropers[order(x@impropers$atm2, x@impropers$atm3), ]
  rownames(x@impropers) <- NULL
  return(x)
}
