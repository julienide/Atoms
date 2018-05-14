#' wrap
#'
#' wrap
#'
#' @name wrap
#' @export
wrap <- function(x, ...)
  UseMethod("wrap")

#' @rdname wrap
#' @export
wrap.Atoms <- function(x, resnumb = x$resnumb, multi = TRUE, ...){
  natm <- natom(x)
  if(is.null(resnumb))
    resnumb <- seq_len(natm)
  if(!is.integer(resnumb))
    stop("'resnumb' must be an integer vector")
  if(length(resnumb) != natm)
    stop("'resnumb' must be of length ", natm)
  resnumb <- as.factor(resnumb)
  if(!is.logical(multi) || length(multi) != 1L)
    stop("'multi' must be a logical vector of length 1")
  if(any(pbc(x))){
    .wrapAtoms(x, resnumb, multi)
  }
  invisible(NULL)
}
