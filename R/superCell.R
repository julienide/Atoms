#' superCell
#' 
#' 
#' @name superCell
#' @export
superCell <- function(x, ...)
  UseMethod("superCell")

#' @rdname superCell
#' @export
superCell.Atoms <- function(x, aInds = 0L, bInds = 0L, cInds = 0L,
                            current = NULL){
  if(is.null(current))
    current <- current(x)
  if(!is.integer(aInds))
    stop("'aInds' must be an integer vector")
  if(!is.integer(bInds))
    stop("'bInds' must be an integer vector")
  if(!is.integer(cInds))
    stop("'cInds' must be an integer vector")
  if(!is.integer(current))
    stop("'current' must be an integer vector")
  
  Obj <- x[, current]
  
  abc <- expand.grid(aInds, bInds, cInds)
  ncell <- nrow(abc)

  frac <- as.matrix(fractional(x))
  natm <- natom(x)

  coords <- rbind(
    rep(frac[, 1L], ncell) + abc[rep(1:ncell, each = natm), 1L],
    rep(frac[, 2L], ncell) + abc[rep(1:ncell, each = natm), 2L],
    rep(frac[, 3L], ncell) + abc[rep(1:ncell, each = natm), 3L])
  coords <- cell(x)%*%coords
  Obj@x <- matrix(coords[1L, ], ncol = 1L)
  Obj@y <- matrix(coords[2L, ], ncol = 1L)
  Obj@z <- matrix(coords[3L, ], ncol = 1L)
  rm(coords)

  Obj@a <- x@a*length(aInds)
  Obj@b <- x@b*length(bInds)
  Obj@c <- x@c*length(cInds)
  Obj@atoms <- Obj@atoms[rep(1:natm, ncell), ]

  duplicateTopo <- function(x, ncell, natm){
    if(nrow(x)){
      x <- x[rep(1:nrow(x), ncell), ] + 
        rep(seq.int(0L, (ncell - 1L)*natm, natm), each = nrow(x))
      rownames(x) <- NULL
    }
    return(x)
  }
  
  Obj@bonds     <- duplicateTopo(Obj@bonds    , ncell, natm)
  Obj@angles    <- duplicateTopo(Obj@angles   , ncell, natm)
  Obj@dihedrals <- duplicateTopo(Obj@dihedrals, ncell, natm)
  Obj@impropers <- duplicateTopo(Obj@impropers, ncell, natm)

  Obj@call <- match.call()
  return(Obj)
}
