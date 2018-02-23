#' neighbors
#' 
#' 
#' @name neighbors
#' @export
neighbors <- function(x, ...)
  UseMethod("neighbors")

#' @rdname neighbors
#' @export
neighbors.Atoms <- function(x, inds, rank = 1L, ...){
  if(missing(inds))
    stop("Please specify 'inds'")
  if(!is.integer(inds))
    stop("'inds' must be an integer vector")
  if(!is.integer(rank) || length(rank) != 1L)
    stop("'rank' must be an integer vector of length 1")
  if(rank < -1L)
    stop("'rank' must be greater or equal to -1")
  if(length(inds)){
    if(rank == 0){
      inds
    } else if(rank == 1){
      ninds <- unique(
        c(
          x@bonds$atm1[x@bonds$atm2 %in% inds],
          x@bonds$atm2[x@bonds$atm1 %in% inds]) )
      ninds[!ninds %in% inds]
    } else if(rank == -1L){
      prev <- 0L
      ninds <- inds
      while(prev != length(ninds)){
        prev <- length(ninds)
        inds <- neighbors(x, inds, rank = 1L)
        ninds <- unique(c(ninds, inds))
      }
      ninds
    } else {
      # For loop not so small here...
      for(n in 1:rank){
        inds <- neighbors(x, inds, rank = 1L)
      }
      inds
    }
  } else {
    integer(0)
  }
}
