#' @rdname Atoms
#' @noRd
#' @export
NCHeader <- function(file, close = TRUE)
{
  con <- nc_open(file, readunlim = FALSE)
  if(close)
    nc_close(con)
  obj <- structure(
    list(
      file = con$filename, con = con,
      natom = con$dim$atom$len, nframe = con$dim$frame$len),
    class = "NCHeader")
  return(obj)
}

#' @rdname Atoms
#' @noRd
#' @export
.NCParser <- function(files, selection = -1L, first = 1L, last = -1L, stride = 1L)
{
  getHdr <- function(file)
  {
    hdr <- NCHeader(file)
    c(natom = hdr$natom, nframe = hdr$nframe,
      cell = !is.null(hdr$con$var$cell_lengths),
      time = !is.null(hdr$con$var$time))
  }
  hdrs <- as.data.frame(t(sapply(files, getHdr)))

  badNatom <- (hdrs$natom != hdrs$natom[1L])
  if(any(badNatom))
    stop("Files with different number of atoms (!= ", hdrs$natom[1L], "): '",
         paste(files[badNatom], collapse = "' "), "'")

  badCell <- hdrs$cell != hdrs$cell[1L]
  if(any(badCell))
    stop("Files with different PBC (!= ", hdrs$cell[1L],"): ",
         paste(files[badCell], collapse = "' "), "'")

  hdrs <- list(
    natom = hdrs$natom[1L],
    nframe = sum(hdrs$nframe),
    cell = hdrs$cell[1L])

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
  
  if(hdrs$cell)
  {
    a <- matrix(NA_real_, nrow = 3L, ncol = nSel)
    b <- matrix(NA_real_, nrow = 3L, ncol = nSel)
    c <- matrix(NA_real_, nrow = 3L, ncol = nSel)
  }
  x <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  y <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  z <- matrix(NA_real_, nrow = length(selection), ncol = nSel)
  
  atomFirst <- min(selection)
  atomCount <- max(selection) - atomFirst + 1L
  # frameFirst <- first
  # frameCount <- last - frameFirst + 1L
  selection <- selection - atomFirst + 1L
  # frameSel <- seq(first, last, stride) - frameFirst + 1L

  sel <- 0L
  pb <- txtProgressBar(0, last, style = 3)
  globalFrame = 0L
  for(file in files)
  {
    hdr <- NCHeader(file, close = FALSE)
    for(frame in 1:hdr$nframe){
      globalFrame <- globalFrame + 1L
      keep <- (globalFrame >= first) && ((globalFrame - first)%%stride == 0)
      err <- try(expr = {
        if(keep){
          sel <- sel + 1L
          if(!is.null(hdr$con$var$cell_lengths)){
            abc <- ncvar_get(hdr$con, "cell_lengths", c(1L, frame), c(3L, 1L))
            abg <- ncvar_get(hdr$con, "cell_angles", c(1L, frame), c(3L, 1L))
            M <- cell2mat(abc[1L], abc[2L], abc[3L], abg[1L], abg[2L], abg[3L])
            a[, sel] <- M[, 1L]
            b[, sel] <- M[, 2L]
            c[, sel] <- M[, 3L]
          }
          start <- c(atomFirst, frame)
          count <- c(1L, atomCount, 1L)
          x[, sel] <- ncvar_get(hdr$con, "coordinates", c(1L, start), count)[selection]
          y[, sel] <- ncvar_get(hdr$con, "coordinates", c(2L, start), count)[selection]
          z[, sel] <- ncvar_get(hdr$con, "coordinates", c(3L, start), count)[selection]
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
    nc_close(hdr$con)
    if(globalFrame == last)
    {
      break
    }
  }
  
  obj <- new("Atoms")
  if(hdrs$cell)
  {
    obj@a <- a ; obj@b <- b ; obj@c <- c ; rm(a, b, c)
    obj@pbc <- rep(TRUE, 3L)
  }
  obj@x <- x ; obj@y <- y ; obj@z <- z ; rm(x, y, z)
  obj@atoms <- data.frame(atmtype = rep("Xx", length(selection)))
  obj@current <- 1L
  
  return(obj)
}
