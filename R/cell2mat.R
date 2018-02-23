#' From Lattice Parameters to Lattice Vectors and Vice Versa
#' 
#' Convert lattice parameters into lattice vectors and vice versa.
#' 
#' 
#' @param a a numeric value. Length of the first lattice vector.
#' @param b a numeric value. Length of the second lattice vector.
#' @param c a numeric value. Length of the third lattice vector.
#' @param alpha a numeric value. Angle between the second and third lattice 
#'   vectors.
#' @param beta a numeric value. Angle between the third and first lattice 
#'   vectors.
#' @param gamma a numeric value. Angle between the first and second lattice 
#'   vectors.
#' @param orient a character string. Orientation of the lattice vectors with 
#'   respect to the Cartesian reference. See details.
#' @param x a 3x3 numeric matrix containing by column the coordinates of the 
#'   lattice vectors.
#'   
#' @return \code{cell2mat} returns a 3x3 numeric matrix containing by column the
#'   coordinates of the lattice vectors.
#'   
#'   \code{mat2cell} returns a data.frame containing by column the lattice
#'   parameters a, b, c, alpha, beta, gamma.
#'   
#' @note Lengths are in Angstrom and angles in degree.
#'   
#' @name cell2mat
#' @export
cell2mat <- function(a = 1L, b = 1L, c = 1L,
                     alpha = 90.0, beta = 90.0, gamma = 90.0,
                     orient = "abxy"){
  if(!is.numeric(a) || length(a) != 1L)
    stop("'a' must be a numeric vector of length 1")
  if(!is.numeric(b) || length(b) != 1L)
    stop("'b' must be a numeric vector of length 1")
  if(!is.numeric(c) || length(c) != 1L)
    stop("'c' must be a numeric vector of length 1")
  if(!is.numeric(alpha) || length(alpha) != 1L)
    stop("'alpha' must be a numeric vector of length 1")
  if(!is.numeric(beta) || length(beta) != 1L)
    stop("'beta' must be a numeric vector of length 1")
  if(!is.numeric(gamma) || length(gamma) != 1L)
    stop("'gamma' must be a numeric vector of length 1")
  if(!is.character(orient) && length(orient) != 1L)
    stop("'orient' must be a character vector of length 1")
  cellVect <- c(substr(orient, 1L, 1L), substr(orient, 2L, 2L))
  cartVect <- c(substr(orient, 3L, 3L), substr(orient, 4L, 4L))
  if(!all(cellVect %in% c("a", "b", "c")) ||
     !all(cartVect %in% c("x", "y", "z")) ||
     cellVect[1L] == cellVect[2L] ||
     cartVect[1L] == cartVect[2L])
    stop("Invalid cell orientation")
  Pabc <- match(cellVect, c("a", "b", "c"))
  Pabc <- c(Pabc, (1:3)[!(1:3 %in% Pabc)])
  Pabc.S <- -sign(Pabc[2] - Pabc[3])
  Pxyz <- match(cartVect, c("x", "y", "z"))
  Pxyz <- c(Pxyz, (1:3)[!(1:3 %in% Pxyz)])
  Pxyz.S <- -sign(Pxyz[2] - Pxyz[3])
  S <- Pabc.S*Pxyz.S
  abc <- c(a, b, c)[Pabc]
  abg <- c(alpha, beta, gamma)[Pabc]*pi/180.0
  V <- prod(abc)*sqrt(1 - sum(cos(abg)^2) + 2*prod(cos(abg)))
  M <- matrix(0.0, nrow = 3L, ncol = 3L)
  M[1L, 1L] <- abc[1L]
  M[1L, 2L] <- abc[2L]*cos(abg[3L])
  M[2L, 2L] <- abc[2L]*sin(abg[3L])
  M[1L, 3L] <- abc[3L]*cos(abg[2L])
  M[2L, 3L] <- abc[3L]*(cos(abg[1L]) - cos(abg[2L])*cos(abg[3L]))/sin(abg[3L])
  M[3L, 3L] <- abc[3L]*sqrt(1 - sum(cos(abg)^2) + 2*prod(cos(abg)))/sin(abg[3L])
    # <- V/(prod(abc[1:2])*sin(abg[3L])) # ! V/prod(abc)
  rownames(M) <- c("x", "y", "z")[Pxyz]
  colnames(M) <- c("a", "b", "c")[Pabc]
  M[3L, ] <- S*M[3L, ]
  M <- M[c("x", "y", "z"), c("a", "b", "c")]
  return(M)
}

#' @name cell2mat
#' @export
mat2cell <- function(x){
  if(!is.matrix(x) || !is.numeric(x) || any(dim(x) != 3L))
    stop("'x' must be a 3x3 numeric matrix")
  obj <- data.frame(
    a = sqrt(x[1L, 1L]^2 + x[2L, 1L]^2 + x[3L, 1L]^2),
    b = sqrt(x[1L, 2L]^2 + x[2L, 2L]^2 + x[3L, 2L]^2),
    c = sqrt(x[1L, 3L]^2 + x[2L, 3L]^2 + x[3L, 3L]^2))
  obj$alpha <- x[1L, 2L]*x[1L, 3L] + x[2L, 2L]*x[2L, 3L] + x[3L, 2L]*x[3L, 3L]
  obj$beta  <- x[1L, 3L]*x[1L, 1L] + x[2L, 3L]*x[2L, 1L] + x[3L, 3L]*x[3L, 1L]
  obj$gamma <- x[1L, 1L]*x[1L, 2L] + x[2L, 1L]*x[2L, 2L] + x[3L, 1L]*x[3L, 2L]
  obj$alpha <- acos(obj$alpha/(obj$b*obj$c))*180/pi
  obj$beta  <- acos(obj$beta /(obj$c*obj$a))*180/pi
  obj$gamma <- acos(obj$gamma/(obj$a*obj$b))*180/pi
  return(obj)
}
