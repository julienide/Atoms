#' Visualize atomic structure
#' 
#' @param \dots rgl.materials
#' 
#' @name view
#' @export
view <- function(x, ...)
  UseMethod("view")

#' @rdname view
#' @export
view.Atoms <- function(
  x, drawing = "lines", color = NULL, radius = NULL, add = FALSE,
  FOV = 0.0, windowRect = c(0.0, 0.0, 800, 600), userMatrix = diag(4),
  pbc.box = NULL, ...){
  # Create representation class for drawing?
  if(is.null(color)){
    if(is.null(x$color)){
      x$color <- atomColor(x)
      message("'atomColor' has been used to guess 'color'")
    }
  } else {
    if(!is.vector(color))
      stop("'color' must be a vector")
    x$color <- color
  }
  if(!all(areColors(x$color)))
    stop("'color' must be a vector of valid R colors")

  if(!add){
    open3d(FOV = FOV, windowRect = windowRect, userMatrix = userMatrix)
    # layout3d(matrix(c(1L, 2L, 1L, 1L), nrow = 2L, ncol = 2L),
    #          heights = windowRect[4L]*c(h <- 0.05, 1 - h),
    #          widths = windowRect[3L]*c(w <- 0.05, 1 - w),
    #          sharedMouse = TRUE)
    rgl.material(...)
  }
  par3d(skipRedraw = TRUE)
  # next3d()
  switch(
    drawing,
    "lines" = {
      if(nrow(bonds(x))){
        inds <- t(bonds(x)[, 1:2])
        segments3d(x$x[inds], x$y[inds], x$z[inds], col = x$color[inds])
      } else {
        warning("No bonds to draw")
        points3d(x$x, x$y, x$z, col = x$color)
      }
    },
    "spheres" = {
      if(is.null(radius)){
        if(is.null(x$radius)){
          x$radius <- rcov(x)
          message("'rcov' has been used to guess 'radius'")
        }
      } else {
        if(!is.numeric(radius))
          stop("'radius' must be a numeric vector")
        x$radius <- radius
      }
      spheres3d(x$x, x$y, x$z, radius = x$radius, col = x$color)
    },
    "points" = {
      points3d(x$x, x$y, x$z, col = x$color)
    },
    {
      stop("Unrocognized drawing style")
    }
  )
  # next3d()
  # # Draw xyz-axes
  if(is.null(pbc.box)){
    if(any(pbc(x)))
      pbc.box <- TRUE
    else
      pbc.box <- FALSE
  }
  if(pbc.box){
    pbcBox(x)
  }
  par3d(skipRedraw = FALSE)
}
