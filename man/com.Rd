% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/com.R
\name{.comAtoms}
\alias{.comAtoms}
\alias{com}
\alias{com.Atoms}
\title{Center Of Mass}
\usage{
.comAtoms(x, mass, factor)

com(x, mass, factor, ...)

\method{com}{Atoms}(x, mass = x$mass, factor = x$resnumb, append = FALSE,
  ...)
}
\arguments{
\item{x}{an R object}

\item{mass}{a numeric vector. atomic masses. If \code{x$mass} is \code{NULL},
function \code{mass} is called to guess masses from atomic symbols.}

\item{factor}{a factor used to group atomic coordinates to compute multiple
centers of mass. If \code{x$resnumb} is \code{NULL} then the center of mass
of the whole structure is computed.}

\item{append}{a logical value. Wheither to append the centers of mass (dummy atoms) to x or
to create a new object.}
}
\value{
an object of class \sQuote{Atoms} containing dummy atoms at each
  centers of mass positions. The masses and residue names of each center of
  mass are also attached to the object.
}
\description{
Calculate the center of mass of multiple residues.
}

