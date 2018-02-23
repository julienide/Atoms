#' join
#' 
#' join
#' 
#' @name join
#' @export
join <- function(x, ...)
  UseMethod("join")

#' @rdname join
#' @export
join.Atoms <- function(x, multi = TRUE, ...){
  if(!is.logical(multi) || length(multi) != 1L)
    stop("'multi' must be a logical vector of length 1")
  if(any(pbc(x))){
    B <- bonds(x)
    B[B$atm1 > B$atm2] <- B[B$atm1 > B$atm2, 2:1]
    B <- B[order(B$atm1, B$atm2), ]
    bonds(x) <- B
    .joinAtoms(x, multi)
  }
  invisible(NULL)
}