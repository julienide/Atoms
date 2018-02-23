#' Lattice Parameters/Vectors
#' 
#' Extract or replace lattice parameters/vectors.
#' 
#' When \code{multi == FALSE}, \code{cell} extracts the lattice vectors of the 
#' current frame. Otherwise, the lattice parameters of all the frames are 
#' extracted.
#' 
#' When \code{multi == FALSE}, \code{cell<-} replace the lattice vectors of the 
#' current frame. Otherwise, the lattice parameters of all the frames are 
#' replaced.
#' 
#' @return When \code{multi = FALSE}, \code{cell} returns a 3x3 numeric matrix 
#'   containing by column the coordinates of the lattice vectors of the current 
#'   frame. Otherwise, a data.frame containing the lattice parameters \code{a}, 
#'   \code{b}, \code{c}, \code{alpha}, \code{beta} and \code{gamma} of each 
#'   frame of the trajectory is returned.
#'   
#'   \code{cell<-} returns an object of the same class as \code{x} with updated 
#'   lattice parameters/vectors.
#'   
#' @param x an R object containing lattice parameters/vectors to 
#'   extracted/replaced.
#' @param multi a logical value. Whether to extract/replace the lattice 
#'   parameters/vectors of all the frames of the trajectory.
#' @param value a 3x3 numeric matrix containing by column the coordinates of the
#'   lattice vectors (when \code{multi == FALSE}) or a data.frame containing by
#'   column the lattice parameters (when \code{multi == TRUE}) to be replaced.
#'   In this case, columns must be named \code{a}, \code{b}, \code{c},
#'   \code{alpha}, \code{beta} or \code{gamma}.
#' 
#' @note Lengths are in Angstrom and angles in degree.
#' 
#' @seealso \code{\link{current}} \code{\link{cell2mat}}
#'   
#' @name cell
#' @export
setGeneric(name = "cell", def = function(x, ...) standardGeneric(f = "cell"))

#' @rdname cell
#' @export
setGeneric(name = "cell<-", def = function(x, ..., value) standardGeneric(f = "cell<-"))

#' @rdname cell
#' @export
setMethod(
  f = "cell",
  signature = "Atoms",
  definition = function(x, multi = FALSE, ...){
    if(!is.logical(multi) || length(multi) != 1L)
      stop("'multi' must be a logical vector of length 1")
    if(multi){
      obj <- data.frame(
        a = sqrt(x@a[1L, ]^2 + x@a[2L, ]^2 + x@a[3L, ]^2),
        b = sqrt(x@b[1L, ]^2 + x@b[2L, ]^2 + x@b[3L, ]^2),
        c = sqrt(x@c[1L, ]^2 + x@c[2L, ]^2 + x@c[3L, ]^2))
      obj$alpha <- x@b[1L, ]*x@c[1L, ] + x@b[2L, ]*x@c[2L, ] + x@b[3L, ]*x@c[3L, ]
      obj$beta  <- x@c[1L, ]*x@a[1L, ] + x@c[2L, ]*x@a[2L, ] + x@c[3L, ]*x@a[3L, ]
      obj$gamma <- x@a[1L, ]*x@b[1L, ] + x@a[2L, ]*x@b[2L, ] + x@a[3L, ]*x@b[3L, ]
      obj$alpha <- acos(obj$alpha/(obj$b*obj$c))*180/pi
      obj$beta  <- acos(obj$beta /(obj$c*obj$a))*180/pi
      obj$gamma <- acos(obj$gamma/(obj$a*obj$b))*180/pi
    } else {
      obj <- cbind(
        x@a[, x@current],
        x@b[, x@current],
        x@c[, x@current])
      colnames(obj) <- c("a", "b", "c")
      rownames(obj) <- c("x", "y", "z")
    }
    return(obj)
  }
)

#' @rdname cell
#' @export
setMethod(
  f = "cell<-",
  signature = "Atoms",
  definition = function(x, multi = FALSE, ..., value){
    if(!is.logical(multi) || length(multi) != 1L)
      stop("'multi' must be a logical vector of length 1")
    if(multi){
      frame <- TRUE
    } else {
      frame <- x@current
    }
    if(length(x@a)){
      if(!is.matrix(value) || !is.numeric(value) || !all(dim(value) == 3L))
        stop("'value' must be a 3x3 numeric matrix")
      x@a[, frame] <- value[, 1L]
      x@b[, frame] <- value[, 2L]
      x@c[, frame] <- value[, 3L]
    } else {
      stop("Empty trajectory")
    }
    return(x)
  }
)
