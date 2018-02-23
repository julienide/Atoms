#' @rdname Atoms
#' @noRd
#' @export
DCDHeader <- function(file, close = TRUE) {
  if(!is.character(file) || length(file) != 1L)
    stop("'file' must be a character vector of length 1")

  con <- file(file, open ="rb")

  # Endianness detection
  end <- .Platform$endian
  FH <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header
  if(FH != 84) # Fortran header should be 84
  {
    file.ori <- summary(con)$description
    close(con)
    end <- ifelse(end == "little", yes = "big", no = "little")
    con <- file(file.ori, "rb")
    FH <- readBin(con, "integer", n = 1L, endian = end)
    if(FH != 84)
    {
      close(con)
      stop("Problem with endian detection")
    }
  }

  # Header
  hdr <- readChar(con, nchars=4) # data => CORD or VELO
  icntrl <- readBin(con, "integer", n = 20L, endian = end)
  names(icntrl) <- c(
    "nfile", "nprevious", "step", "nstep", rep("", 3L),
    "ndegf", "nfrozen", "delta", "crystal", "fourDims", rep("", 7L), "charmmVersion")
  FE <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header

  if(icntrl["charmmVersion"])
  {
    # re-read X-PLOR delta as a double
    cur.pos <- seek(con, where = 44, origin = "start")
    icntrl["delta"]  <- readBin(con, "double", 1, endian = end)
    seek(con, where = cur.pos, origin = "start")
  }

  # Titles
  FH <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header
  ntitle <- readBin(con, "integer", 1, endian = end)
  cur.pos <- seek(con, where = NA)  # 100
  title <- rep(NA_character_, ntitle)
  for(n in seq_along(title))
    title[n] <- try(readChar(con, 80), silent = TRUE)
  if(any(inherits(title, "try-error")))
    stop("Invalid title")
  FE <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header

  # Number of atoms
  FH <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header
  natom <- readBin(con, "integer", 1, endian = end)
  FE <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header

  # Free (movable) atom indexes
  if(icntrl["nfrozen"])
  {
    FH <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header
    free <- readBin(con, "integer", (natom - icntrl["nfrozen"]), endian = end)
    FE <- readBin(con, "integer", n = 1L, endian = end) # FORTRAN header
  } else {
    free = integer(0)
  }

  # Close file if requested
  if(close)
    close(con)

  obj <- structure(
    list(
      file = file, con = con, endianness = end,
      hdr = hdr, icntrl = icntrl, title = title,
      free = free, natom = as.integer(natom),
      nframe = as.integer(unname(icntrl["nfile"]))),
    class = "DCDHeader")
  return(obj)
}

#' #' @rdname Atoms
#' #' @noRd
#' #' @export
#' show.DCDHeader <- function(object)
#' {
#'   show(object[names(object) != "con"])
#' }
#' 
#' #' @rdname Atoms
#' #' @noRd
#' #' @export
#' print.DCDHeader <- function(x, ...)
#' {
#'   print(x[names(x) != "con"])
#' }

#' @rdname Atoms
#' @noRd
#' @export
.DCDParser <- function(files, selection = -1L, first = 1L, last = -1L, stride = 1L)
{
  readCell <- function(hdr)
  {
    FH <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
    cell <- readBin(hdr$con, "numeric", n = (FH/8), endian = hdr$endianness, size = 8)
    FE <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
    return(cell[c(1L, 3L, 6L, 5L, 4L, 2L)])
  }

  readXYZ <- function(hdr)
  {
    FH <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
    x <- readBin(hdr$con, "numeric", n=(FH/4), endian = hdr$endianness, size=4)
    FE <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
    return(x)
  }

  skipBlock <- function(hdr)
  {
    FH <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
    seek(hdr$con, where = FH, origin = "current")
    FE <- readBin(hdr$con, "integer", n = 1L, endian = hdr$endianness) # FORTRAN header
  }
  
  getHdr <- function(file)
  {
    hdr <- DCDHeader(file)
    c(hdr$icntrl["crystal"], natom = hdr$natom, nframe = hdr$nframe)
  }
  hdrs <- as.data.frame(t(sapply(files, getHdr)))

  badNatom <- (hdrs$natom != hdrs$natom[1L])
  if(any(badNatom))
    stop("Files with different number of atoms (!= ", hdrs$natom[1L], "): '",
         paste(files[badNatom], collapse = "' "), "'")
  
  badCrystal <- hdrs$crystal != hdrs$crystal[1L]
  if(any(badCrystal))
    stop("Files with different PBC (!= ", hdrs$crystal[1L],"): ",
         paste(files[badCrystal], collapse = "' "), "'")
  
  hdrs <- list(
    natom = hdrs$natom[1L],
    nframe = sum(hdrs$nframe),
    crystal = hdrs$crystal[1L])
  
  if(length(selection) == 1L && selection == -1L)
    selection <- seq_len(hdrs$natom)
  
  if(last == -1L)
  {
    last <- hdrs$nframe
  }
  else
  {
    if(last > hdrs$nframe)
      stop("'last' is out of range (last > ", hdrs$nframe, ")")
  }
    

  nSel <- ceiling((last - first + 1)/stride)

  if(hdrs$crystal)
  {
    a <- matrix(NA_real_, nrow = 3L, ncol = nSel)
    b <- matrix(NA_real_, nrow = 3L, ncol = nSel)
    c <- matrix(NA_real_, nrow = 3L, ncol = nSel)
  }
  x <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  y <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  z <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  
  sel <- 0L
  pb <- txtProgressBar(1, last, style = 3)
  globalFrame = 0L
  for(file in files)
  {
    hdr <- DCDHeader(file, close = FALSE)
    for(frame in 1:hdr$nframe){
      globalFrame <- globalFrame + 1L
      keep <- (globalFrame >= first) && ((globalFrame - first)%%stride == 0)
      err <- try(expr = {
        if(keep){
          sel <- sel + 1L
          if(hdr$icntrl["crystal"])
          {
            cell <- readCell(hdr)
            sinGamma <- sqrt(1 - cell[6L]^2)
            a[, sel] <- cell[1L]*c(1.0, 0.0, 0.0)
            b[, sel] <- cell[2L]*c(cell[6L], sinGamma, 0.0)
            c[, sel] <- cell[3L]*c(cell[5L], (cell[4L] - cell[5L]*cell[6L])/sinGamma,
                                   sqrt(1 - sum(cell[4:6]^2) + 2*prod(cell[4:6]))/sinGamma)
          }
          x[, sel] <- readXYZ(hdr)[selection]
          y[, sel] <- readXYZ(hdr)[selection]
          z[, sel] <- readXYZ(hdr)[selection]
        } else {
          if(hdr$icntrl["crystal"])
          {
            skipBlock(hdr)
          }
          skipBlock(hdr)
          skipBlock(hdr)
          skipBlock(hdr)
        }
        # CHARMM files may contain an extra block
        if (hdr$icntrl["charmmVersion"] && hdr$icntrl["fourDims"]) {
          skipBlock(hdr)
        }
      }, silent = TRUE)
      if(inherits(err, "try-error")){
        stop("Erro while reading file", file," at frame ", frame, "/", hdr$nframe)
      }
      setTxtProgressBar(pb, globalFrame)
      if(globalFrame == last)
      {
        break
      }
    }
    close(hdr$con)
    if(globalFrame == last)
    {
      break
    }
  }
  
  obj <- new("Atoms")
  if(hdrs$crystal)
  {
    obj@a <- a ; obj@b <- b ; obj@c <- c ; rm(a, b, c)
    obj@pbc <- rep(TRUE, 3L)
  }
  obj@x <- x ; obj@y <- y ; obj@z <- z ; rm(x, y, z)
  obj@atoms <- data.frame(atmtype = rep("Xx", length(selection)))
  obj@current <- 1L

  return(obj)
}
