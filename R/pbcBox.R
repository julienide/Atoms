#' draw PBC box
#' 
#' @name pbcBox
#' @export
pbcBox <- function(x, ...)
  UseMethod("pbcBox")

#' @rdname pbcBox
#' @export
pbcBox.matrix <- function(x, ...){
  ## TODO: Draw pbc box for 2D and 1D systems
  if(!is.numeric(x))
    stop("'x' must be a numeric matrix")
  if(nrow(x) != 3L || ncol(x) != 3L)
    stop("'x' must be a 3x3 matrix")
  vertices <- cbind(rep(0.0, 3L), x[, 1L], x[, 1L] + x[, 2L], x[, 2L])
  vertices <- cbind(vertices, vertices + x[, 3L])
  indices <- cbind(
    rbind(1:4, c(2:4, 1L)),
    rbind(5:8, c(6:8, 5L)),
    rbind(1:4, 5:8))
  col <- rep(1L, length(indices))
  col[c(1:2, 7:8, 17:18)] <- rep(1L + 1:3, each = 2L)
  segments3d(t(vertices)[indices, ], col = col, ...)
}

#' @rdname pbcBox
#' @export
pbcBox.Atoms <- function(x, ...){
  pbcBox(cell(x), ...)
}