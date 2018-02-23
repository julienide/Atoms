#' selectAtoms
#' 
#' 
#' @name selectAtoms
#' @export
selectAtoms <- function(x, expr)
  UseMethod("selectAtoms")

#' @rdname selectAtoms
#' @export
selectAtoms.Atoms <- function(x, expr)
{
  x <- c(
    list(
       x = x$x ,  y = x$y ,  z = x$z,
      vx = x$vx, vy = x$vy, vz = x$vz,
      fx = x$fx, fy = x$fy, fz = x$fz),
    as.list(x@atoms))
  x <- x[!is.null(x)]
  eval(substitute(expr), x, enclos = parent.frame())
}