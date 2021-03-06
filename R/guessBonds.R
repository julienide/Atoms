#' #' Guess Bonds
#' #' 
#' #' @name guessBonds
#' #' @export
#' guessBonds <- function(x, ...)
#'   UseMethod("guessBonds")
#' 
#' #' #' @rdname guessBonds
#' #' #' @export
#' #' guessBonds.Atoms <- function(x, radius = NULL, safety = 1.2){
#' #'   if(!is.numeric(safety) || length(safety) != 1L)
#' #'     stop("'safety' must be a numeric vector of length 1")
#' #' 
#' #'   if(is.null(radius)){
#' #'     if(is.null(x$radius)){
#' #'       radius <- rcov(x)
#' #'     } else {
#' #'       radius <- x$radius
#' #'     }
#' #'   } else {
#' #'     if(!is.numeric(radius))
#' #'       stop("'radius' must be a numeric vector")
#' #'     if(length(radius) == 1L)
#' #'       radius <- rep(radius, natom(x))
#' #'     if(length(radius) != natom(x))
#' #'       stop("'radius' must be of length ", natom(x))
#' #'   }
#' #' 
#' #'   guessBondsLocal <- function(x, cell, cellInv, pbc){
#' #'     if(nrow(x) > 1L){
#' #'       nat <- nrow(x)
#' #'       dx <- outer(x$x, x$x, "-")
#' #'       dy <- outer(x$y, x$y, "-")
#' #'       dz <- outer(x$z, x$z, "-")
#' #' 
#' #'       h <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
#' #'       k <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
#' #'       l <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
#' #' 
#' #'       # Not possible to find bond going throught the cell with this approach
#' #'       # as bond are computed within small boxes only and not inbetween them
#' #'       if(any(pbc)){
#' #'         da <- dx*cellInv[1L, 1L] + dy*cellInv[1L, 2L] + dz*cellInv[1L, 3L]
#' #'         db <- dx*cellInv[2L, 1L] + dy*cellInv[2L, 2L] + dz*cellInv[2L, 3L]
#' #'         dc <- dx*cellInv[3L, 1L] + dy*cellInv[3L, 2L] + dz*cellInv[3L, 3L]
#' #'         # print(da[da > 0.5 | da < -0.5])
#' #'         if(pbc[1L]){
#' #'           h <- floor(da + 0.5)
#' #'           da <- da - h
#' #'         }
#' #'         if(pbc[2L]){
#' #'           k <- floor(db + 0.5)
#' #'           db <- db - k
#' #'         }
#' #'         if(pbc[3L]){
#' #'           l <- floor(dc + 0.5)
#' #'           dc <- dc - l
#' #'         }
#' #'         dx <- da*cell[1L, 1L] + db*cell[1L, 2L] + dc*cell[1L, 3L]
#' #'         dy <- da*cell[2L, 1L] + db*cell[2L, 2L] + dc*cell[2L, 3L]
#' #'         dz <- da*cell[3L, 1L] + db*cell[3L, 2L] + dc*cell[3L, 3L]
#' #'       }
#' #' 
#' #'       ds <- dx*dx + dy*dy + dz*dz
#' #'       dbond <- outer(x$r, x$r, "+")
#' #'       dbond <- dbond*dbond
#' #'       M <- lower.tri(ds) & (0.1 < ds) & (ds < dbond)
#' #' 
#' #'       inds <- matrix(x$atmnumb, nrow = nat, ncol = nat)
#' #'       data.frame(
#' #'         atm1 = t(inds)[M], atm2 = inds[M],
#' #'         h = h[M], k = k[M], l = l[M])
#' #'     } else {
#' #'       data.frame(
#' #'         atm1 = integer(0), atm2 = integer(0),
#' #'         h = integer(0), k = integer(0), l = integer(0))
#' #'     }
#' #'   }
#' #' 
#' #'   guessBonds <- function(shift, step, x, cell, cellInv, pbc){
#' #'     xrange <- range(x$x, na.rm = TRUE) + shift[1L]
#' #'     yrange <- range(x$y, na.rm = TRUE) + shift[2L]
#' #'     zrange <- range(x$z, na.rm = TRUE) + shift[3L]
#' #' 
#' #'     xcut <- cut(x$x, seq(xrange[1L], xrange[2L] + step, step), include.lowest = TRUE)
#' #'     ycut <- cut(x$y, seq(yrange[1L], yrange[2L] + step, step), include.lowest = TRUE)
#' #'     zcut <- cut(x$z, seq(zrange[1L], zrange[2L] + step, step), include.lowest = TRUE)
#' #' 
#' #'     seq(xrange[1L])
#' #'     
#' #'     do.call(
#' #'       rbind,
#' #'       by(x, list(xcut, ycut, zcut),
#' #'          guessBondsLocal, cell, cellInv, pbc, simplify = FALSE) )
#' #'   }
#' #' 
#' #'   cell <- cell(x)
#' #'   cellInv <- solve(cell)
#' #'   pbc <- pbc(x)
#' #'   x <- data.frame(x = x$x, y = x$y, z = x$z,
#' #'                   r = radius*safety, atmnumb = 1:natom(x))
#' #' 
#' #'   step <- 10.0
#' #'   shift <- expand.grid(0:1,0:1,0:1)*step/2L
#' #' 
#' #'   B <- do.call(
#' #'     rbind,
#' #'     apply(shift, 1L, guessBonds, step, x, cell, cellInv, pbc) )
#' #' 
#' #'   B <- unique(B)
#' #'   B[order(B$atm1, B$atm2), ]
#' #' }