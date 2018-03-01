#' Center Of Mass
#'
#' Calculate the center of mass of multiple residues.
#'
#' @param x an R object
#' @param mass a numeric vector. atomic masses. If \code{x$mass} is \code{NULL},
#'   function \code{mass} is called to guess masses from atomic symbols.
#' @param factor a factor used to group atomic coordinates to compute multiple
#'   centers of mass. If \code{x$resnumb} is \code{NULL} then the center of mass
#'   of the whole structure is computed.
#' @param append a logical value. Wheither to append the centers of mass (dummy atoms) to x or
#'   to create a new object.
#'
#' @return an object of class \sQuote{Atoms} containing dummy atoms at each
#'   centers of mass positions. The masses and residue names of each center of
#'   mass are also attached to the object.
#'
#' @name com
#' @export
com <- function(x, mass, factor, ...)
  UseMethod("com")

#' @rdname com
#' @export
com.Atoms <- function(x, mass = x$mass, factor = x$resnumb, append = FALSE, ...){
  if(!length(x@x))
    stop("'x' doesn't contain coordinates")

  if(is.null(mass)){
    mass <- mass(x)
    message("'mass' has been used to guess 'mass'")
  }
  if(!is.numeric(mass))
    stop("'mass' must be a numeric vector")
  if(length(mass) != natom(x))
    stop("'mass' must be of length ", natom(x))

  if(is.null(factor))
    factor <- rep(1L, natom(x))
  if(!is.factor(factor))
    factor <- as.factor(factor)
  if(length(factor) != natom(x))
    stop("'factor' must be of length ", natom(x))

  if(!is.logical(append) || length(append) != 1L)
    stop("'append' must be a logical vector of length 1")

  com <- .comAtoms(x, mass, factor)
  if(append){
    x@x <- rbind(x@x, com@x)
    x@y <- rbind(x@y, com@y)
    x@z <- rbind(x@z, com@z)
    natm <- natom(x)
    ncom <- natom(com)
    x@atoms <- x@atoms[c(1:natm, rep(NA_integer_, ncom)), ]
    rownames(x@atoms) <- NULL
    if(!is.null(x$atmtype)){
      x@atoms$atmtype[natm + 1:ncom] <- "Xx"
    }
    if(!is.null(x$atmname)){
      x@atoms$atmname[natm + 1:ncom] <- levels(factor)
    }
    if(!is.null(x$mass)){
      x@atoms$mass[natm + 1:ncom] <- com$mass
    }
    return(x)
  } else {
    return(com)
  }

  # Obj <- new("Atoms")
  # Obj@pbc <- x@pbc
  # Obj@current <- 1L
  # Obj@a <- x@a
  # Obj@b <- x@b
  # Obj@c <- x@c
  # Obj@call <- match.call()
  # Obj@atoms <- data.frame(
  #   atmname = "Xx",
  #   atmtype = "Xx",
  #   resname = levels(factor),
  #   mass = tapply(mass, factor, sum),
  #   stringsAsFactors = FALSE)
  #
  # mass <- mass/unsplit(sapply(split(mass, factor), sum), factor)
  # Obj@x <- x@x*mass
  # Obj@y <- x@y*mass
  # Obj@z <- x@z*mass
  #
  # nframe <- nframe(x)
  # factor <- as.factor(outer(factor, seq(nframe), paste))
  # Obj@x <- matrix(tapply(Obj@x, factor, sum), ncol = nframe)
  # Obj@y <- matrix(tapply(Obj@y, factor, sum), ncol = nframe)
  # Obj@z <- matrix(tapply(Obj@z, factor, sum), ncol = nframe)
  #
  # return(Obj)
}

