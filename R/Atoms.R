#' Create an object of class \sQuote{Atoms}
#' 
#' This function read topology and trajectory files to build an object of class 
#' \sQuote{Atoms}.
#' 
#' To build an object of class \sQuote{Atoms} a topology and trajectory file are
#' required. A batch of files can also also be provided for the trajectory. For 
#' the sake of simplicity, when \code{topo.file}/\code{trj.files} is missing, 
#' the function attempts to read the topology/trajectory from 
#' \code{trj.files}/\code{topo.file}.
#' 
#' By default, topology and trajectory file formats are guessed from file
#' extensions. Otherwise, they can be specified with the \code{topo.fmt} and 
#' \code{trj.fmt} arguments. Currently, the following file formats can be read: 
#' XYZ, PDB (ENT), RES, DCD, NC (NETCDF, NCDF) and LAMMPS.
#' 
#' The reading can be restricted to specific atoms and frames by using the 
#' \code{selection}, \code{first}, \code{last} and \code{stride} arguments.
#' 
#' @param topo.file a character vector. Topology file name. See details.
#' @param trj.files a character vector. Trajectory file name(s). See details.
#' @param selection an integer vector. Indexes of atoms to be read.
#' @param first an integer. First frame to be read.
#' @param last an integer. last frame to be read.
#' @param stride an integer. Step used for reading the trajectory frames.
#' @param topo.fmt a character string. Format used for reading topology file. 
#'   See details
#' @param trj.fmt a character string. Format used for reading trajectory file. 
#'   See details
#'
#' @include guessFMT.R
#' @importFrom Rcpp evalCpp
#' @useDynLib Atoms
#'   
#' @name Atoms
#' @export
Atoms <- function(topo.file, trj.files, selection = -1L,
                  first = 1L, last = -1L, stride = 1L,
                  topo.fmt = NULL, trj.fmt = NULL)
{
  if(missing(topo.file) || !length(topo.file)){
    if(missing(trj.files) || !length(trj.files)){
      stop("Please specify at least one topology or trajectory file")
    } else {
      topo.file <- trj[1L]
    }
  } else {
    if(missing(trj.files) || !length(trj.files)){
      trj.files <- topo.file
      topo.file <- topo.file[1L]
    } else {
      if(length(topo.file) != 1L){
        warning("Only the first topology file will be read")
        topo.file <- topo.file[1L]
      }
    }
  }
  
  if(!is.character(topo.file))
    stop("'topo.file' must be a character vector")
  topo.file <- path.expand(topo.file)
  if(!file.exists(topo.file))
    stop("File '", topo.file, "' is missing")

  if(!is.character(trj.files))
    stop("'trj.files' must be a character vector")
  trj.files <- path.expand(trj.files)
  missingFiles <- !file.exists(trj.files)
  if(any(missingFiles))
    stop("Missing files: '", paste(trj.files[missingFiles], sep = "', "), "'")
  
  if(!is.null(topo.fmt)){
    if(!is.character(topo.fmt) || length(topo.fmt) != 1L)
      stop("'topo.fmt' must be a character vector of length 1")
  }
  if(!is.null(trj.fmt)){
    if(!is.character(trj.fmt) || length(trj.fmt) != 1L)
      stop("'trj.fmt' must be a character vector of length 1")
  }
  
  if(!is.integer(selection) && !is.logical(selection))
    stop("'selection' must be logical or integer vector")
  
  if(anyNA(selection))
    stop("'NA' value are not permitted in 'selection'")
  
  if(is.integer(selection))
  {
    if(any(selection < 1L))
    {
      if(length(selection) != 1L)
      {
        stop("'selection' must contain atom indices (>= 1L)")
      } else {
        if(selection != -1L)
          stop("'selection' must be -1 (all atoms) or contain atom indices (!= 0)")
      }
    }
    if(is.unsorted(selection))
    {
      selection <- sort(selection)
      warning("'selection' has been sorted")
    }
  }
  if(is.logical(selection))
  {
    if(length(selection) == 1L)
    {
      if(selection)
      {
        selection <- -1L
      }
      else
      {
        stop("Please selecte at least one atom")
      }
    }
    else
    {
      selection <- which(selection)
    }
  }
  if(length(selection) == 0)
    stop("Empty atomic selection")
  
  if(!is.integer(first) || length(first) != 1L)
    stop("'first' must be an integer vector of length 1")
  if(!is.integer(last) || length(last) != 1L)
    stop("'last' must be an integer vector of length 1")
  if(!is.integer(stride) || length(stride) != 1L)
    stop("'stride' must be an integer vector of length 1")
  
  if(stride < 1L)
    stop("'stride' must be greater than zero")
  if(first < 1L)
    stop("'first' must be greater than zero")
  if(last < 1L && last != -1L)
    stop("'last' must equal to -1L (all frames) or greater than zero")
  if(last != -1L && first > last)
    stop("'first' is greater than 'last'")

  if(is.null(topo.fmt))
  {
    topo.fmt <- guessFMT(topo.file)
    message(topo.fmt, " topology file format has been guessed")
  }
  topo.parser <- paste0(".", topo.fmt, "Parser")
  if(!exists(topo.parser))
    stop("No parser implemented for '", topo.fmt, "' file format")
  
  if(is.null(trj.fmt))
  {
    trj.fmt <- guessFMT(trj.files)
    message(trj.fmt, " trajectory file format has been guessed")
  }
  trj.parser <- paste0(".", trj.fmt, "Parser")
  if(!exists(trj.parser))
    stop("No parser implemented for '", trj.fmt, "' file format")
  
  trj <- do.call(
    trj.parser,
    list(trj.files, selection = selection,
         first = first, last = last, stride = stride))
  if(!topo.fmt %in% c("DCD", "NC")){
    topology(trj) <- do.call(
      topo.parser,
      list(topo.file, selection = selection, last = 1L))
  }
  
  trj@call <- match.call()
  validObject(trj)
  
  return(trj)
}

