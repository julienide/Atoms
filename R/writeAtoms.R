#' writeAtoms
#' 
#' writeAtoms
#' 
#' @name writeAtoms
#' @export
writeAtoms <- function(x, file, fmt = NULL, multi = FALSE){
  if(!is.character(file) || length(file) != 1L)
    stop("'file' must be a character string of length 1")
  file <- path.expand(file)
  if(is.null(fmt)){
    fmt <- guessFMT(file)
    message(fmt, " file format has been guessed")
  }
  if(!is.character(fmt) || length(fmt) != 1L)
    stop("'fmt' must be a character vector of length 1")
  writer <- paste0(".", fmt, "Writer")
  if(!exists(writer))
    stop("No writer implemented for '", fmt, "' file format")
  if(fmt == "RES"){
    multi <- TRUE
  }
  if(multi){
    file <- strsplit(file, "\\.")[[1L]]
    ext <- file[length(file)]
    file <- paste(file[-length(file)], collapse = ".")
    fmt <- paste0("%s_%0", nchar(nframe(x)), "i.%s")
    file <- sprintf(fmt = fmt, file, seq_len(nframe(x)), ext)
  }
  # else {
  #   file <- rep(file, nframe(x))
  # }
  x@bonds <- x@bonds[order(x@bonds$atm1, x@bonds$atm2), ]
  do.call(writer, list(x, file))
}