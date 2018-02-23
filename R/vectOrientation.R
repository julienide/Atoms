#' Vectors Orientation
#' 
#' Calculate the angles between two sets of vectors (Distance vectors or plan 
#' normals).
#' 
#' When \code{V1} and \code{V2} are data.frame, they must contain columns named 
#' \code{atm1}, \code{atm2} and optionally \code{atm3}. They must all be integer
#' vectors containing atom indexes defining the vectors used for the analysis. 
#' When only \code{atm1} and \code{atm2} are specified, distances vectors are 
#' used. When \code{atm3} is also specified, plan normals are used instead. In
#' this latter case, plan normals are calculated as the cross product between
#' two normalized distance vectors d12 (defined by \code{atm1} and \code{atm2})
#' and d13 (defined by \code{atm1} and \code{atm3}).
#' 
#' When \code{V2} is missing, the angles between V1 and itself are calculated 
#' (lower-triangle of the (V1, V1) angles matrix). When \code{V2} is a numeric 
#' vector, the angles between V1 and this single reference vector are 
#' calculated. When both \code{V1} and \code{V2} are data.frame, the angles 
#' between V1 and V2 are calculated (a (N1,N2) matrix for each frame is 
#' calculated).
#' 
#' @return A 2D or 3D numeric array. The last dimension span over the frames of 
#'   the trajectory.
#'   
#' @param V1 a data.frame containing pairs/triplets of atom indexes defining a 
#'   set of distance/normal vectors.
#' @param V2 Same as V1 or a numeric vector of length 3 indicating a single 
#'   reference vector for the calculation of the angles. See details.
#' @param pairwise a logical value. When both V1 and V2 are data.frame, 
#'   indicates wheither to compute angles pairwise or not.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @name vectOrientation
#' @export
vectOrientation <- function(x, ...)
  UseMethod("vectOrientation")

#' @rdname vectOrientation
#' @export
vectOrientation.Atoms <- function(x, V1, V2, pairwise = TRUE, ...){
  if(!is.data.frame(V1))
    stop("'V1' must be a data.frame")
  if(is.null(V1$atm1))
    stop("'V1' must contain a 'atm1' column")
  if(is.null(V1$atm2))
    stop("'V1' must contain a 'atm2' column")
  atm3 <- V1$atm3
  V1 <- data.frame(
    atm1 = checkSel(V1$atm1, natom(x)),
    atm2 = checkSel(V1$atm2, natom(x)))
  if(!is.null(atm3)){
    V1$atm3 <- checkSel(atm3, natom(x)) ; rm(atm3)
  }

  if(missing(V2)){
    .vectOrientation(x, V1)
  } else {
    if(is.numeric(V2)){
      if(length(V2) != 3L)
        stop("'V2' must be of length 3")
      .vectOrientationAxis(x, V1, V2)
    } else {
      if(!is.data.frame(V2))
        stop("'V2' must be a data.frame")
      if(is.null(V2$atm1))
        stop("'V2' must contain a 'atm1' column")
      if(is.null(V2$atm2))
        stop("'V2' must contain a 'atm2' column")
      atm3 <- V2$atm3
      V2 <- data.frame(
        atm1 = checkSel(V2$atm1, natom(x)),
        atm2 = checkSel(V2$atm2, natom(x)))
      if(!is.null(atm3)){
        V2$atm3 <- checkSel(atm3, natom(x)) ; rm(atm3)
      }
      if(pairwise){
        if(nrow(V1) != nrow(V2))
          stop("'V1' and 'V2' must have the same number of rows for pairwise calculation")
        .vectOrientationPairwise(x, V1, V2)
      } else {
        .vectOrientationNotPairwise(x, V1, V2)
      }
    }
  }
}