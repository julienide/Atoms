#' Volume and Density of a Simulation Cell
#' 
#' Calculate the volume or density of a simulation cell.
#' 
#' @param x an R object
#' @param multi a logical value. Whether to calculate the volume/density of the 
#'   simulation cell for all the frames of the trajectory.
#' @param mass a numeric vector. atomic masses. If \code{x$mass} is \code{NULL},
#'   function \code{mass} is called to guess masses from atomic symbols.
#'   
#' @return \code{cellVolume} returns a numeric vector containing the volume of
#'   the simulation cell.
#'   
#'   \code{cellDensity} returns a numeric vector containing the density of the
#'   simulation cell.
#'   
#' @note The volume is in Angstrom^3 cube and the density in g/cm^2.
#'   
#' @name cellProperties
#' @export
setGeneric(name = "cellVolume", def = function(x, ...) standardGeneric(f = "cellVolume"))

setMethod(
  f = "cellVolume",
  signature = "matrix",
  definition = function(x, ...){
    dims <- dim(x)
    if(!is.numeric(x) || length(dims) != 2L || any(dims != 3L))
      stop("'x' must be a 3x3 numeric matrix")
    V <- 
      x[1L, 1L]*(x[2L, 2L]*x[3L, 3L] - x[3L, 2L]*x[2L, 3L]) +
      x[1L, 2L]*(x[3L, 2L]*x[1L, 3L] - x[1L, 2L]*x[3L, 3L]) +
      x[1L, 3L]*(x[1L, 2L]*x[2L, 3L] - x[2L, 2L]*x[1L, 3L])
    return(V)
  }
)

#' @rdname cellProperties
#' @export
setMethod(
  f = "cellVolume",
  signature = "Atoms",
  definition = function(x, multi = FALSE, ...){
    if(!is.logical(multi) || length(multi) != 1L)
      stop("'multi' must be a logical vector of length 1")
    if(multi){
      frame <- TRUE
    } else {
      frame <- x@current
    }
    V <- 
      x@a[1L, frame]*(
        x@b[2L, frame]*x@c[3L, frame] - x@b[3L, frame]*x@c[2L, frame]) +
      x@a[2L, frame]*(
        x@b[3L, frame]*x@c[1L, frame] - x@b[1L, frame]*x@c[3L, frame]) +
      x@a[3L, frame]*(
        x@b[1L, frame]*x@c[2L, frame] - x@b[2L, frame]*x@c[1L, frame])
    return(V)
  }
)

#' @name cellProperties
#' @export
setGeneric(name = "cellDensity", def = function(x, ...) standardGeneric(f = "cellDensity"))

#' @rdname cellProperties
#' @export
setMethod(
  f = "cellDensity",
  signature = "Atoms",
  definition = function(x, mass = x$mass, multi = FALSE, ...){
    if(is.null(mass)){
      mass <- mass(x)
      message("'mass' has been used to guess 'mass'")
    }
    if(!is.numeric(mass))
      stop("'mass' must be a numeric vector")
    if(length(mass) != natom(x))
      stop("'mass' must be of length ", natom(x))
    D <- sum(mass, na.rm = TRUE)/cellVolume(x, multi = multi)
    D <- D/(1E-24*6.0221412927E+23)
    return(D)
  }
)

