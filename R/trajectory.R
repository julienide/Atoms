#' trajectory
#' 
#' Get or set the lattice parameters, atomic poistions, velocities and forces of
#' an object.
#' 
#' @param x an R object from/to which to get/set the trajectory
#' @param value an R object containing the trajectory used for setting \code{x}
#'   
#' @name trajectory
#' @export
trajectory <- function(x)
  UseMethod("trajectory")

trajectory.Atoms <- function(x){
  x@atoms <- data.frame(atmtype = rep("Xx", natom(x)))
  x@bonds <- x@bonds[integer(0), ]
  x@angles <- x@angles[integer(0), ]
  x@dihedrals <- x@dihedrals[integer(0), ]
  x@impropers <- x@impropers[integer(0), ]
  x@call <- match.call()
  return(x)
}

"trajectory<-" <- function(x, value)
  UseMethod("trajectory<-")

"trajectory<-.Atoms" <- function(x, value){
  # Generally value is smaller than x (contain the trajectory)
  # Then it is less memory consuming to set the topology of 'value' and return it
  # than setting the trajectory of 'x' and return it.
  if(natom(x) != natom(value))
    stop("'x' and 'value' contain a different number of atoms")
  topology(value) <- topology(x)
  current(value) <- current(x)
  value@call <- match.call()
  return(value)
}