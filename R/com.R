#' Center Of Mass
#' 
#' Calculate the center of mass multiple residues.
#' 
#' @param x an R object
#' @param mass a numeric vector. atomic masses. If \code{x$mass} is \code{NULL},
#'   function \code{mass} is called to guess masses from atomic symbols.
#' @param factor a factor used to group atomic coordinates to compute multiple 
#'   centers of mass. If \code{x$resnumb} is \code{NULL} then the center of mass
#'   of the whole structure is computed.
#'   
#' @return an object of class \sQuote{Atoms} containing dummy atoms at each 
#'   centers of mass positions. The masses and residue names of each center of
#'   mass are also attached to the object.
#'   
#' @name com
#' @export
com <- function(x, mass, factor)
  UseMethod("com")

#' @rdname com
#' @export
com.Atoms <- function(x, mass = x$mass, factor = x$resnumb){
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
  
  .comAtoms(x, mass, factor)
  
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

