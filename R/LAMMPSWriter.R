#' writeLAMMPS
#' 
#' 
#' @name writeLAMMPS
#' @export
writeLAMMPS <- function(x, file, type)
  UseMethod("writeLAMMPS")

#' @rdname writeLAMMPS
#' @export
writeLAMMPS.Atoms <- function(x, file, type = "full"){
  LAMMPSDataTypes <- c(
    "angle", "atomic", "body", "bond",
    "charge", "dipole", "dpd", "electron",
    "ellipsoid", "full", "line", "meso",
    "molecular", "peri", "smd", "sphere",
    "template", "tri", "wavepacket", "hybrid")
  
  if(!is.character(file) || length(file) != 1L)
    stop("'file' must be a cahracter vector of length 1")
  if(!is.character(type) || length(type) != 1L)
    stop("'type' must be a cahracter vector of length 1")
  if(!type %in% LAMMPSDataTypes)
    stop("'type' must be one of: ", paste(LAMMPSDataTypes, collapse = ", "))
  
  if(is.null(x$atmtype))
    stop("Missing atom-type ('x$atmtype' is null)")
  atmtype <- as.integer(as.factor(x$atmtype))

  atoms <- switch(
    type,
    angle = {stop("Not implemented yet")},
    atomic = {
      # atom-ID atom-type x y z
      sprintf(
        fmt = "%6i %6i %10.5f %10.5f %10.5f # %s",
        1:natom(x), atmtype, x$x, x$y, x$z, x$atmtype)
    },
    body = {stop("Not implemented yet")},
    bond = {stop("Not implemented yet")},
    charge = {
      if(is.null(x$charge))
        stop("Missing charges ('x$charge' is null)")
      # atom-ID atom-type q x y z
      sprintf(
        fmt = "%6i %6i %10.5f %10.5f %10.5f %10.5f # %s",
        1:natom(x), atmtype, x$charge, x$x, x$y, x$z, x$atmtype)
    },
    dipole = {stop("Not implemented yet")},
    dpd = {stop("Not implemented yet")},
    electron = {stop("Not implemented yet")},
    ellipsoid = {stop("Not implemented yet")},
    full = {
      if(is.null(x$resnumb))
        stop("Missing molecule-ID ('x$resnumb' is null)")
      if(is.null(x$charge))
        stop("Missing charges ('x$charge' is null)")
      # atom-ID molecule-ID atom-type q x y z
      sprintf(
        fmt = "%6i %6i %6i %10.5f %10.5f %10.5f %10.5f # %s",
        1:natom(x), x$resnumb, atmtype, x$charge, x$x, x$y, x$z, x$atmtype)
    },
    line = {stop("Not implemented yet")},
    meso = {stop("Not implemented yet")},
    molecular = {stop("Not implemented yet")},
    peri = {stop("Not implemented yet")},
    smd = {stop("Not implemented yet")},
    sphere = {stop("Not implemented yet")},
    template = {stop("Not implemented yet")},
    tri = {stop("Not implemented yet")},
    wavepacket = {stop("Not implemented yet")},
    hybrid = {stop("Not implemented yet")}
  )
  
  if(is.null(x$atmtype))
    x$atmtype <- x$atmname
  
  fmt <- max(nchar(x$atmtype))
  fmt <- paste0("%-", fmt, "s")
  
  types <- topoDict(x, unique = FALSE)
  types$bonds <- as.factor(
    sprintf(
      fmt = paste(rep(fmt, 2L), collapse = " "),
      types$bonds$atm1,
      types$bonds$atm2))
  types$angles <- as.factor(
    sprintf(
      fmt = paste(rep(fmt, 3L), collapse = " "),
      types$angles$atm1,
      types$angles$atm2,
      types$angles$atm3))
  types$dihedrals <- as.factor(
    sprintf(
      fmt = paste(rep(fmt, 4L), collapse = " "),
      types$dihedrals$atm1,
      types$dihedrals$atm2,
      types$dihedrals$atm3,
      types$dihedrals$atm4))
  types$impropers <- as.factor(
    sprintf(
      fmt = paste(rep(fmt, 4L), collapse = " "),
      types$impropers$atm1,
      types$impropers$atm2,
      types$impropers$atm3,
      types$impropers$atm4))
  
  btypes <- as.integer(types$bonds)
  atypes <- as.integer(types$angles)
  dtypes <- as.integer(types$dihedrals)
  itypes <- as.integer(types$impropers)
  
  types$bonds <- unique(types$bonds)
  types$angles <- unique(types$angles)
  types$dihedrals <- unique(types$dihedrals)
  types$impropers <- unique(types$impropers)
  
  types$bonds <- structure(as.integer(types$bonds), names = as.character(types$bonds))
  types$angles <- structure(as.integer(types$angles), names = as.character(types$angles))
  types$dihedrals <- structure(as.integer(types$dihedrals), names = as.character(types$dihedrals))
  types$impropers <- structure(as.integer(types$impropers), names = as.character(types$impropers))
  
  cell <- round(cell(x), digits = 5L)
  
  lines <- c("LAMMPS data file", "")
  lines <- c(lines, sprintf(fmt = "%6i atoms", natom(x)))
  if(nrow(bonds(x)))
    lines <- c(lines, sprintf(fmt = "%6i bonds", nrow(bonds(x))))
  if(nrow(angles(x)))
    lines <- c(lines, sprintf(fmt = "%6i angles", nrow(angles(x))))
  if(nrow(dihedrals(x)))
    lines <- c(lines, sprintf(fmt = "%6i dihedrals", nrow(dihedrals(x))))
  if(nrow(impropers(x)))
    lines <- c(lines, sprintf(fmt = "%6i impropers", nrow(impropers(x))))
  lines <- c(lines, "")
  cat(lines, file = file, sep = "\n")

  lines <- character(0)
  if(length(types$bonds))
    lines <- c(lines, sprintf(fmt = "%6i bond types", length(types$bonds)))
  if(length(types$angles))
    lines <- c(lines, sprintf(fmt = "%6i angle types", length(types$angles)))
  if(length(types$dihedrals))
    lines <- c(lines, sprintf(fmt = "%6i dihedral types", length(types$dihedrals)))
  if(length(types$impropers))
    lines <- c(lines, sprintf(fmt = "%6i improper types", length(types$impropers)))
  lines <- c(lines, "")
  cat(lines, file = file, sep = "\n", append = TRUE)

  lines <- c(
    sprintf(fmt = "%10.5f %10.5f xlo xhi", 0.0, cell[1L, 1L]),
    sprintf(fmt = "%10.5f %10.5f ylo yhi", 0.0, cell[2L, 2L]),
    sprintf(fmt = "%10.5f %10.5f zlo zhi", 0.0, cell[3L, 3L]))
  if(any(cell[upper.tri(cell)] != 0.0)){
    lines <- c(
      lines,
      sprintf(
        fmt = "%10.5f %10.5f %10.5f xy xz yz",
        cell[1L, 2L], cell[1L, 3L], cell[2L, 3L]))
  }
  cat(lines, "", file = file, sep = "\n", append = TRUE)
  
  if(length(types$bonds)){
    cat(
      "Bond Coeffs",
      sprintf(fmt = "%6i    0.0000    0.0000 # %s",
              types$bonds, names(types$bonds)),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(length(types$angles)){
    cat(
      "Angle Coeffs",
      sprintf(fmt = "%6i    0.0000    0.0000 # %s",
              types$angles, names(types$angles)),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(length(types$dihedrals)){
    cat(
      "Dihedrals Coeffs",
      sprintf(fmt = "%6i    0.0000    0.0000 # %s",
              types$dihedrals, names(types$dihedrals)),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(length(types$impropers)){
    cat(
      "Impropers Coeffs",
      sprintf(fmt = "%6i    0.0000    0.0000 # %s",
              types$impropers, names(types$impropers)),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  
  cat("Atoms",
      atoms,
      "",
      file = file, sep = "\n", append = TRUE)
  
  if(nrow(bonds(x))){
    cat(
      "Bonds",
      sprintf(fmt = "%6i %6i %6i %6i",
              1:nrow(bonds(x)), btypes, bonds(x)$atm1, bonds(x)$atm2),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(nrow(angles(x))){
    cat(
      "Angles",
      sprintf(fmt = "%6i %6i %6i %6i %6i",
              1:nrow(angles(x)), atypes,
              angles(x)$atm1, angles(x)$atm2, angles(x)$atm3),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(nrow(dihedrals(x))){
    cat(
      "Dihedrals",
      sprintf(fmt = "%6i %6i %6i %6i %6i %6i",
              1:nrow(dihedrals(x)), dtypes,
              dihedrals(x)$atm1, dihedrals(x)$atm2,
              dihedrals(x)$atm3, dihedrals(x)$atm4),
      "",
      file = file, sep = "\n", append = TRUE)
  }
  if(nrow(impropers(x))){
    cat(
      "Impropers",
      sprintf(fmt = "%6i %6i %6i %6i %6i %6i",
              1:nrow(impropers(x)), itypes,
              impropers(x)$atm1, impropers(x)$atm2,
              impropers(x)$atm3, impropers(x)$atm4),
      "",
      file = file, sep = "\n", append = TRUE)
  }
}
