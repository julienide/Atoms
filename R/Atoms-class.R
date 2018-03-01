#' Atoms Class
#'
#' An S4 class to store the topology and trajectory of an atomic system.
#'
#' Atomic positions, velocities and forces as well as lattice parameters of
#' multiple frames of a molecular dynamics trajectory can be stored into an
#' \sQuote{Atoms} object. The topology of the system (bonds, bond angles,
#' dihedrals and impropers angles) can also be stored. Any type of atomic
#' properties can easily be attached to the object.
#'
#' @section \sQuote{\code{$}} and \sQuote{\code{$<-}} methods: Atomic properties
#'   stored inside slot \code{atoms} can be get or set using the dollard syntax.
#'   Atomic positions (\code{x}, \code{y} and \code{z}), velocities (\code{vx},
#'   \code{vy} and \code{vz}) and forces (\code{fx}, \code{fy} and \code{fz}) of
#'   the current frame can also be get or set.
#'
#' @section \sQuote{\code{[}} and \sQuote{\code{[<-}} methods: Extract or
#'   replace a subset of an object of class \sQuote{Atoms}. The first/second
#'   index span over the atoms/frames. \code{NA} values can be used for atom
#'   selection to introduce gap positions. \code{NA} values are not permited for
#'   frame selection.
#'
#' @slot atoms a \code{data.frame} containg atomic properties.
#' @slot bonds a \code{data.frame} containing the folling columns: \describe{
#'   \item{atm1,atm2}{integer vectors. Atom indexes forming bonds.} }
#' @slot angles a \code{data.frame} containing the folling columns: \describe{
#'   \item{atm1,atm2,atm3}{integer vectors. Atom indexes forming bond angles.} }
#' @slot dihedrals a \code{data.frame} containing the folling columns:
#'   \describe{ \item{atm1,atm2,atm3, atm4}{integer vectors. Atom indexes
#'   forming dihedral angles.} }
#' @slot impropers a \code{data.frame} containing the folling columns:
#'   \describe{ \item{atm1,atm2,atm3, atm4}{integer vectors. Atom indexes
#'   forming improper angles.} }
#' @slot current an integer vector of length 1. Index of the current frame.
#' @slot pbc a logical vector of length 3. Periodic boundary conditions along
#'   the three directions of the cell.
#' @slot a,b,c numeric matrices. Cartesian coordinates (rows) of lattice vectors
#'   for each frame (columns).
#' @slot x,y,z numeric matrices. Cartesian coordinates (rows) of atomic position
#'   vectors for each frame (columns).
#' @slot vx,vy,vz numeric matrices. Cartesian coordinates (rows) of velocity
#'   vectors for each frame (columns).
#' @slot fx,fy,fz numeric matrices. Cartesian coordinates (rows) of force
#'   vectors for each frame (columns).
#' @slot call the call used to create the object.
#'
#' @keywords classes
#'
#' @name Atoms-class
#' @export
setClass(
  Class = "Atoms",
  slots = list(
    atoms = "data.frame",
    bonds = "data.frame",
    angles = "data.frame",
    dihedrals = "data.frame",
    impropers = "data.frame",
    current = "integer",
    pbc = "logical",
    a = "matrix", b = "matrix", c = "matrix",
    x = "matrix", y = "matrix", z = "matrix",
    vx = "matrix", vy = "matrix", vz = "matrix",
    fx = "matrix", fy = "matrix", fz = "matrix",
    call = "call" ),
  prototype = list(
    bonds = data.frame(
      atm1 = integer(0),
      atm2 = integer(0)),
    angles = data.frame(
      atm1 = integer(0),
      atm2 = integer(0),
      atm3 = integer(0)),
    dihedrals = data.frame(
      atm1 = integer(0),
      atm2 = integer(0),
      atm3 = integer(0),
      atm4 = integer(0)),
    impropers = data.frame(
      atm1 = integer(0),
      atm2 = integer(0),
      atm3 = integer(0),
      atm4 = integer(0)),
    current = 0L,
    pbc = rep(FALSE, 3L),
    a = matrix(numeric(0), nrow = 3L),
    b = matrix(numeric(0), nrow = 3L),
    c = matrix(numeric(0), nrow = 3L),
    x = matrix(numeric(0), nrow = 0L, ncol = 0L),
    y = matrix(numeric(0), nrow = 0L, ncol = 0L),
    z = matrix(numeric(0), nrow = 0L, ncol = 0L),
    vx = matrix(numeric(0), nrow = 0L, ncol = 0L),
    vy = matrix(numeric(0), nrow = 0L, ncol = 0L),
    vz = matrix(numeric(0), nrow = 0L, ncol = 0L),
    fx = matrix(numeric(0), nrow = 0L, ncol = 0L),
    fy = matrix(numeric(0), nrow = 0L, ncol = 0L),
    fz = matrix(numeric(0), nrow = 0L, ncol = 0L) ),
  validity = function(object){
    checkColumn <- function(slotname, colname, object){
      if(is.null(slot(object, slotname)[[colname]])){
        paste0("Slot '", slotname, "': Missing '", colname, "' column")
      } else if(!is.integer(slot(object, slotname)[[colname]])){
        paste0("Slot '", slotname, "': Column '", colname, "' must be integer")
      } else{
        TRUE
      }
    }
    msg <- checkColumn("bonds", "atm1", object)
    msg <- checkColumn("bonds", "atm2", object)
    msg <- checkColumn("angles", "atm1", object)
    msg <- checkColumn("angles", "atm2", object)
    msg <- checkColumn("angles", "atm3", object)
    msg <- checkColumn("dihedrals", "atm1", object)
    msg <- checkColumn("dihedrals", "atm2", object)
    msg <- checkColumn("dihedrals", "atm3", object)
    msg <- checkColumn("dihedrals", "atm4", object)
    msg <- checkColumn("impropers", "atm1", object)
    msg <- checkColumn("impropers", "atm2", object)
    msg <- checkColumn("impropers", "atm3", object)
    msg <- checkColumn("impropers", "atm4", object)
    if(is.character(msg))
      return(msg)
    natom <- nrow(object@atoms)
    inds <- seq(natom)
    if(!all(slot(object, "bonds")$atm1 %in% inds) ||
       !all(slot(object, "bonds")$atm2 %in% inds) )
      return("Slot 'bonds' contains indices out of atom range [1, ", natom,"]")
    if(!all(slot(object, "angles")$atm1 %in% inds) ||
       !all(slot(object, "angles")$atm2 %in% inds) ||
       !all(slot(object, "angles")$atm3 %in% inds) )
      return("Slot 'angles' contains indices out of atom range [1, ", natom,"]")
    if(!all(slot(object, "dihedrals")$atm1 %in% inds) ||
       !all(slot(object, "dihedrals")$atm2 %in% inds) ||
       !all(slot(object, "dihedrals")$atm3 %in% inds) ||
       !all(slot(object, "dihedrals")$atm4 %in% inds) )
      return("Slot 'dihedrals' contains indices out of atom range [1, ", natom,"]")
    if(!all(slot(object, "impropers")$atm1 %in% inds) ||
       !all(slot(object, "impropers")$atm2 %in% inds) ||
       !all(slot(object, "impropers")$atm3 %in% inds) ||
       !all(slot(object, "impropers")$atm4 %in% inds) )
      return("Slot 'impropers' contains indices out of atom range [1, ", natom,"]")
    if(length(slot(object, "current")) != 1L)
      return("Slot 'current' must be of length 1")
    if(length(slot(object, "pbc")) != 3L)
      return("Slot 'pbc' must be of length 3")
    if(!is.numeric(slot(object, "a")))
      return("Slot 'a' must be a numeric matrix")
    if(!is.numeric(slot(object, "b")))
      return("Slot 'b' must be a numeric matrix")
    if(!is.numeric(slot(object, "c")))
      return("Slot 'c' must be a numeric matrix")
    if(!is.numeric(slot(object, "x")))
      return("Slot 'x' must be a numeric matrix")
    if(!is.numeric(slot(object, "y")))
      return("Slot 'y' must be a numeric matrix")
    if(!is.numeric(slot(object, "z")))
      return("Slot 'z' must be a numeric matrix")
    if(!is.numeric(slot(object, "vx")))
      return("Slot 'vx' must be a numeric matrix")
    if(!is.numeric(slot(object, "vy")))
      return("Slot 'vy' must be a numeric matrix")
    if(!is.numeric(slot(object, "vz")))
      return("Slot 'vz' must be a numeric matrix")
    if(!is.numeric(slot(object, "fx")))
      return("Slot 'fx' must be a numeric matrix")
    if(!is.numeric(slot(object, "fy")))
      return("Slot 'fy' must be a numeric matrix")
    if(!is.numeric(slot(object, "fz")))
      return("Slot 'fz' must be a numeric matrix")
    if(nrow(slot(object, "a")) != 3L)
      return("Slot 'a' must a 3xN matrix")
    if(nrow(slot(object, "b")) != 3L)
      return("Slot 'b' must a 3xN matrix")
    if(nrow(slot(object, "c")) != 3L)
      return("Slot 'c' must a 3xN matrix")
    dims <- data.frame(
      a = dim(slot(object, "a")),
      b = dim(slot(object, "b")),
      c = dim(slot(object, "c")),
      x = dim(slot(object, "x")),
      y = dim(slot(object, "y")),
      z = dim(slot(object, "z")),
      vx = dim(slot(object, "x")),
      vy = dim(slot(object, "y")),
      vz = dim(slot(object, "z")),
      fx = dim(slot(object, "x")),
      fy = dim(slot(object, "y")),
      fz = dim(slot(object, "z")))
    if(!all(dims$a == dims$b) || !all(dims$a == dims$c))
      return("Slots 'a', 'b' and 'c' must have the same dimensions")
    if(!all(dims$x == dims$y) || !all(dims$x == dims$z))
      return("Slots 'x', 'y' and 'z' must have the same dimensions")
    if(!all(dims$vx == dims$vy) || !all(dims$vx == dims$vz))
      return("Slots 'vx', 'vy' and 'vz' must have the same dimensions")
    if(!all(dims$fx == dims$fy) || !all(dims$fx == dims$fz))
      return("Slots 'fx', 'fy' and 'fz' must have the same dimensions")
    if(any(dims[2L, ] != 0L)){
      if(dims$a[2L] == 0L || dims$b[2L] == 0L || dims$b[2L] == 0L)
        return(
          paste0(
            "Slots 'a', 'b', 'c' must be defined when ",
            "'x', 'y', 'z', ",
            "'vx', 'vy', 'vz', ",
            "'vx', 'vy' or 'vz' is defined"))
      if(any(dims[2L, ] != 0L & dims[2L, ] != dims$a[2L]))
        return(
          paste0(
            "Slots ",
            "'a', 'b', 'c', ",
            "'x', 'y', 'z', ",
            "'vx', 'vy', 'vz', ",
            "'vx', 'vy' and 'vz' ",
            "must have the same number of columns"))
    }
    nframe <- dims$a[2L]
    nrows <- unlist(dims[1L, -(1:3)])
    if(any(nrows != 0L)){
      nrows <- nrows[nrows != 0L]
      if(any(nrows != nrows[1L]))
        return(
          paste0(
            "Slots ",
            "'x', 'y', 'z', ",
            "'vx', 'vy', 'vz', ",
            "'vx', 'vy' and 'vz' ",
            "must have the same number of rows"))
    }
    natom <- nrows[1L]
    current <- slot(object, "current")
    if(current == 0L){
      if(natom != 0L || nframe != 0L)
        return("Slot 'current' can be equal to 0 only when the object is empty")
    } else {
      if(current < 1L || current > nframe)
        return(paste0("Slot 'current' is out of range [1, ", nframe,"]"))
    }
    if(natom != 0L && natom != nrow(object@atoms))
      return(paste0("Slots 'atoms' must contain the same number of rows as ",
                    "'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy' and 'fz'"))
    TRUE
  }
)

#' @rdname Atoms-class
#' @export
setMethod(
  f = "$",
  signature = "Atoms",
  definition = function(x, name){
    getable <- c("x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")
    if(name %in% getable){
      if(length(slot(x, name)))
        slot(x, name)[, x@current]
      else
        NULL
    } else {
      x@atoms[[name]]
    }
  }
)

## TODO: Replacement of x, y or z values by NA values must be done consitantly for each coordinate
#' @rdname Atoms-class
#' @export
setMethod(
  f = "$<-",
  signature = "Atoms",
  definition = function(x, name, value){
    setable <- c("x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")
    if(name %in% setable){
      if(!is.numeric(value))
        stop("'value' must be a numeric vector")
      if(length(value) != 1L && length(value) != natom(x))
        stop("replacement has ", length(value), " rows, data has ", natom(x))
      if(length(slot(x, name))){
        slot(x, name)[, x@current] <- value
      } else {
        slot(x, name) <- matrix(value, nrow = natom(x), ncol = 1L)
        def <- gsub("v||f", "", name)
        undef <- grep(def, c("x", "y", "z"), invert = TRUE, value = TRUE)
        for(u in undef){
          uname <- gsub(def, u, name)
          slot(x, uname) <- matrix(0.0, nrow = natom(x), ncol = 1L)
        }
        slot(x, "current") <- 1L
      }

      # if(length(slot(x, name))){
      #   if(length(value) != 1L && length(value) != natom(x))
      #     stop("replacement has ", length(value), " rows, data has ", natom(x))
      #   slot(x, name)[, x@current] <- value
      # } else {
      #   stop("Can not replace '", name, "' which is not defined")
      # }
    } else {
      x@atoms[[name]] <- value
    }
    return(x)
  }
)

#' @rdname Atoms-class
#' @export
setMethod(
  f = "[",
  signature = "Atoms",
  definition = function(x, i, j, ..., drop = FALSE){
    if(missing(i))
      i <- TRUE
    if(missing(j))
      j <- TRUE

    if(!is.logical(i) && !is.numeric(i))
      stop("Selection must be done with logical or integer vectors")
    if(!is.logical(j) && !is.numeric(j))
      stop("Selection must be done with logical or integer vectors")
    if(anyNA(j))
      stop("NA values are not permitted for frame selection")

    natom <- natom(x)
    if(is.logical(i))
      i <- seq_len(natom)[i]
    else
      i <- as.integer(i)
    i[i > natom] <- NA # Impose same behavious as a data.frame for atom indexing.

    nframe <- nframe(x)
    if(is.logical(j))
      j <- seq_len(nframe)[j]
    else
      j <- as.integer(j)
    if(any(j <= 0L | j > nframe) || anyNA(j))
      stop("Selection is out of frame range [1, ", nframe, "]")

    x@a <- x@a[, j, drop = FALSE]
    x@b <- x@b[, j, drop = FALSE]
    x@c <- x@c[, j, drop = FALSE]
    if(length(x@x)){
      x@x <- x@x[i, j, drop = FALSE]
      x@y <- x@y[i, j, drop = FALSE]
      x@z <- x@z[i, j, drop = FALSE]
    }
    if(length(x@vx)){
      x@vx <- x@vx[i, j, drop = FALSE]
      x@vy <- x@vy[i, j, drop = FALSE]
      x@vz <- x@vz[i, j, drop = FALSE]
    }
    if(length(x@fx)){
      x@fx <- x@fx[i, j, drop = FALSE]
      x@fy <- x@fy[i, j, drop = FALSE]
      x@fz <- x@fz[i, j, drop = FALSE]
    }
    if(x@current %in% j){
      print(is(x@current))
      print(is(as.integer(min(j))))
      # Shift index of current frame
      x@current <- x@current - as.integer(min(j)) + 1L
    } else {
      # warning("Current frame has been removed. Current frame has been set to 1")
      x@current <- 1L
    }

    x@atoms <- x@atoms[i, , drop = FALSE]
    i <- na.omit(i)
    x@bonds <- x@bonds[
      (x@bonds$atm1 %in% i) &
      (x@bonds$atm2 %in% i), ]
    x@angles <- x@angles[
      (x@angles$atm1 %in% i) &
      (x@angles$atm2 %in% i) &
      (x@angles$atm3 %in% i), ]
    x@dihedrals <- x@dihedrals[
      (x@dihedrals$atm1 %in% i) &
      (x@dihedrals$atm2 %in% i) &
      (x@dihedrals$atm3 %in% i) &
      (x@dihedrals$atm4 %in% i), ]
    x@impropers <- x@impropers[
      (x@impropers$atm1 %in% i) &
      (x@impropers$atm2 %in% i) &
      (x@impropers$atm3 %in% i) &
      (x@impropers$atm4 %in% i), ]
    indexes <- rownames(x@atoms)
    x@bonds$atm1 <- match(x@bonds$atm1, indexes)
    x@bonds$atm2 <- match(x@bonds$atm2, indexes)
    x@angles$atm1 <- match(x@angles$atm1, indexes)
    x@angles$atm2 <- match(x@angles$atm2, indexes)
    x@angles$atm3 <- match(x@angles$atm3, indexes)
    x@dihedrals$atm1 <- match(x@dihedrals$atm1, indexes)
    x@dihedrals$atm2 <- match(x@dihedrals$atm2, indexes)
    x@dihedrals$atm3 <- match(x@dihedrals$atm3, indexes)
    x@dihedrals$atm4 <- match(x@dihedrals$atm4, indexes)
    x@impropers$atm1 <- match(x@impropers$atm1, indexes)
    x@impropers$atm2 <- match(x@impropers$atm2, indexes)
    x@impropers$atm3 <- match(x@impropers$atm3, indexes)
    x@impropers$atm4 <- match(x@impropers$atm4, indexes)
    rownames(x@atoms) <- NULL
    rownames(x@bonds) <- NULL
    rownames(x@angles) <- NULL
    rownames(x@dihedrals) <- NULL
    rownames(x@impropers) <- NULL

    return(x)
  }
)

#' @rdname Atoms-class
#' @export
.DollarNames.Atoms <- function(x, pattern = ""){
  names <- character(0)
  if(length(x@x))
    names <- c(names, "x", "y", "z")
  if(length(x@vx))
    names <- c(names, "vx", "vy", "vz")
  if(length(x@fx))
    names <- c(names, "fx", "fy", "fz")
  names <- c(names(x@atoms), names)
  grep(pattern, names, value = TRUE)
}

