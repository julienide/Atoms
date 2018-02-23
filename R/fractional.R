#' Fraction Coordinates
#' 
#' Extract fractional coordinates of a given frame of a trajectory.
#' 
#' @param x an R object.
#' @param frame an integer. The frame for which to extract the fractional
#'   coordinates.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @return A data.frame containing by columns the xyz fractional coordinates.
#'   
#' @name fractional
#' @export
fractional <- function(x, ...)
  UseMethod("fractional")

#' @rdname fractional
#' @export
fractional.Atoms <- function(x, frame = current(x), ...)
{
  if(!length(x@x))
    stop("'x' doesn't contain coordinates")
  if(!any(pbc(x)))
    stop("'x' doesn't have PBC")
  current(x) <- frame
  cell <- solve(cell(x))
  data.frame(
    x =
      x@x[, frame]*cell[1L, 1L] +
      x@y[, frame]*cell[1L, 2L] +
      x@z[, frame]*cell[1L, 3L],
    y =
      x@x[, frame]*cell[2L, 1L] +
      x@y[, frame]*cell[2L, 2L] +
      x@z[, frame]*cell[2L, 3L],
    z =
      x@x[, frame]*cell[3L, 1L] +
      x@y[, frame]*cell[3L, 2L] +
      x@z[, frame]*cell[3L, 3L])
}


