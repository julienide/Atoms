#' rdf
#' 
#' rdf
#' 
#' @name rdf
#' @include checkSel.R
#' @export
rdf <- function(x, ...)
  UseMethod("rdf")

#' @rdname rdf
#' @export
rdf.Atoms <- function(x, sel1 = NULL, sel2 = NULL, resnumb = x$resnumb,
                      cutoff = 20.0, interval = 0.02, type = "xyz"){
  if(!length(x@x))
    stop("'x' doesn't contain coordinates")
  
  if(is.null(sel1))
    sel1 <- 1:natom(x)
  if(is.null(sel2))
    sel2 <- sel1
  natm <- natom(x)
  sel1 <- checkSel(sel1, natm)
  sel2 <- checkSel(sel2, natm)
  if(is.null(resnumb))
    resnumb <- rep(1L, natm)
  if(!is.integer(resnumb))
    stop("'resnumb' must be an integer vector")
  if(length(resnumb) != natm)
    stop("'resnumb' must be of length ", natm)
  if(!is.numeric(cutoff) || length(cutoff) != 1L)
    stop("'cutoff' must be a numeric vector of length 1")
  if(!is.numeric(interval) || length(interval) != 1L)
    stop("'interval' must be a numeric vector of length 1")
  if(cutoff <= 0.0)
    stop("'cutoff' must be greater than 0")
  if(interval <= 0.0)
    stop("'interval' must be greater than 0")
  type <- match.arg(type, c("xyz", "xy", "yz", "zx", "x", "y", "z"))
  # Calculate Viso and Vloc
  # For xy, yz and zx assum a, b and c to be along x, y, z!!
  # Add a criteria to check the alignment of a, b and c with respect to x, y and z?
  breaks <- seq(0.0, cutoff, interval)
  type <- c(grepl("x", type), grepl("y", type),grepl("z", type))

  Obj <- .rdfAtoms(x, sel1, sel2, resnumb, cutoff, interval, type)
  
  return(Obj)
}
