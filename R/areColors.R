#' Are the elements of a vector valid R colors?
#'
#' Check if the elements of a vector are valid R colors.
#'
#' @return a logical vector indicating if the elements of \code{x} are valid
#'   colors.
#'
#' @examples
#' areColors(c(1, "1", NA, "black", "#000000"))
#'
#' @keywords utilities
#'
#' @seealso \code{\link[grDevices]{colors}}
#'
#' @param x an atomic vector to be tested.
#' @export
areColors <- function(x) {
  if(!is.vector(x))
    stop("'x' must be a vector")
  sapply(x, function(X) {
    tryCatch(
      if(is.na(X))
        NA
      else
        is.matrix(col2rgb(X)),
      error = function(e) FALSE)
  })
}
