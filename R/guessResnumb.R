#' guessResnumb
#' 
#' 
#' @name guessResnumb
#' @export
guessResnumb <- function(x)
  UseMethod("guessResnumb")

#' @rdname guessResnumb
#' @export
guessResnumb.Atoms <- function(x){
  resnumb <- rep(NA_integer_, natom(x))
  N <-  0L
  while(anyNA(resnumb)) {
    N <- N + 1L
    ind <- which(is.na(resnumb))[1L]
    # ind <- which.max( is.na(resnumb))
    # ind <- which.min(!is.na(resnumb))
    resnumb[neighbors(x, ind, rank = -1L)] <- N
  }
  return(resnumb)
}
