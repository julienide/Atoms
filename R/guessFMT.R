guessFMT <- function(file){
  fmt <- tail(strsplit(file, "\\.")[[1L]], 1L)
  fmt <- toupper(fmt)
  if(fmt == "ENT")
    fmt <- "PDB"
  if(fmt %in% c("NETCDF", "NCDF"))
    fmt <- "NC"
  return(fmt)
}