checkSel <- function(x, natom){
  xname <-match.call()$x
  if(!is.integer(x) && !is.logical(x))
    stop("'", xname,"' must be an integer or a logical vector")
  if(is.logical(x)){
    if(length(x) != natom)
      stop("'", xname, "' must be of length ", natom)
    x <- which(x)
  } else if(is.integer(x)){
    if(any(x < 1) || any(x > natom))
      stop("'", xname, "' contains atom indexes out of range")
  } else {
    stop("'", xname, "' must be a logical or integer vector")
  }
  return(x)
}