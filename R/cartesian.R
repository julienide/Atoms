#' Cartesian Coordinates
#'
#' Extract Cartesian coordinates of a given frame of a trajectory.
#'
#' @param x an R object.
#' @param frame an integer. The frame for which to extract the Cartesian
#'   coordinates.
#' @param \dots further arguments passed to or from other methods.
#'
#' @return A data.frame containing by columns the xyz Cartesian coordinates.
#'
#' @name cartesian
#' @export
cartesian <- function(x, ...)
  UseMethod("cartesian")

#' @rdname cartesian
#' @export
cartesian.Atoms <- function(x, frame = current(x), ...)
{
  if(!length(x@x))
    stop("'x' doesn't contain coordinates")
  current(x) <- frame
  data.frame(x = x@x[, frame], y = x@y[, frame], z = x@z[, frame])
}


