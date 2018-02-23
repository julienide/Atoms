#' topology
#' 
#' 
#' @name topology
#' @export
topology <- function(x)
  UseMethod("topology")

#' @rdname topology
#' @export
topology.Atoms <- function(x){
  # Return an object of class Atoms containing only the topology
  x@pbc <- rep(FALSE, 3L)
  x@current <- 0L
  x@a <- matrix(numeric(0), nrow = 3L, ncol = 0L)
  x@b <- matrix(numeric(0), nrow = 3L, ncol = 0L)
  x@c <- matrix(numeric(0), nrow = 3L, ncol = 0L)
  x@x <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  x@y <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  x@z <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  return(x)
}

#' @rdname topology
#' @export
"topology<-" <- function(x, value)
  UseMethod("topology<-")

#' @rdname topology
#' @export
"topology<-.Atoms" <- function(x, value){
  if(!inherits(value, "Atoms"))
    stop("'value' must be an object of class 'Atoms'")
  x@atoms <- value@atoms
  x@bonds <- value@bonds
  x@angles <- value@angles
  x@dihedrals <- value@dihedrals
  x@impropers <- value@impropers
  return(x)
}
