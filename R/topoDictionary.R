#' Topology Dictionary
#'
#' Topology Dictionary
#'
#'
#' @name topoDictionary
#' @export
topoDict <- function(x, atmprop, unique)
  UseMethod("topoDict")

#' @rdname topoDictionary
#' @export
topoDict.Atoms <- function(x, atmprop = "atmtype", unique = TRUE){
  if(is.null(x@atoms[[atmprop]]))
    stop("Undefined atomic property: '", atmprop, "'")
  else
    atmprop <- x@atoms[[atmprop]]
  dict <- list(
    bonds = data.frame(
      atm1 = atmprop[x@bonds$atm1],
      atm2 = atmprop[x@bonds$atm2],
      stringsAsFactors = FALSE),
    angles = data.frame(
      atm1 = atmprop[x@angles$atm1],
      atm2 = atmprop[x@angles$atm2],
      atm3 = atmprop[x@angles$atm3],
      stringsAsFactors = FALSE),
    dihedrals = data.frame(
      atm1 = atmprop[x@dihedrals$atm1],
      atm2 = atmprop[x@dihedrals$atm2],
      atm3 = atmprop[x@dihedrals$atm3],
      atm4 = atmprop[x@dihedrals$atm4],
      stringsAsFactors = FALSE),
    impropers = data.frame(
      atm1 = atmprop[x@impropers$atm1],
      atm2 = atmprop[x@impropers$atm2],
      atm3 = atmprop[x@impropers$atm3],
      atm4 = atmprop[x@impropers$atm4],
      stringsAsFactors = FALSE) )
  if(unique){
    dict$bonds[dict$bonds$atm1 > dict$bonds$atm2, ] <-
      dict$bonds[dict$bonds$atm1 > dict$bonds$atm2, 2:1]
    dict$angles[dict$angles$atm1 > dict$angles$atm3, ] <-
      dict$angles[dict$angles$atm1 > dict$angles$atm3, 3:1]
    dict$dihedrals[dict$dihedrals$atm2 > dict$dihedrals$atm3, ] <-
      dict$dihedrals[dict$dihedrals$atm2 > dict$dihedrals$atm3, 4:1]
    dict$impropers[dict$impropers$atm2 > dict$impropers$atm3, ] <-
      dict$impropers[dict$impropers$atm2 > dict$impropers$atm3, 4:1]
    dict$bonds <- unique(dict$bonds)
    dict$angles <- unique(dict$angles)
    dict$dihedrals <- unique(dict$dihedrals)
    dict$impropers <- unique(dict$impropers)
    row.names(dict$bonds) <- NULL
    row.names(dict$angles) <- NULL
    row.names(dict$dihedrals) <- NULL
    row.names(dict$impropers) <- NULL
  }
  return(dict)
}
