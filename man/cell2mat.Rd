% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell2mat.R
\name{cell2mat}
\alias{cell2mat}
\alias{mat2cell}
\title{From Lattice Parameters to Lattice Vectors and Vice Versa}
\usage{
cell2mat(a = 1L, b = 1L, c = 1L, alpha = 90, beta = 90, gamma = 90,
  orient = "abxy")

mat2cell(x)
}
\arguments{
\item{a}{a numeric value. Length of the first lattice vector.}

\item{b}{a numeric value. Length of the second lattice vector.}

\item{c}{a numeric value. Length of the third lattice vector.}

\item{alpha}{a numeric value. Angle between the second and third lattice 
vectors.}

\item{beta}{a numeric value. Angle between the third and first lattice 
vectors.}

\item{gamma}{a numeric value. Angle between the first and second lattice 
vectors.}

\item{orient}{a character string. Orientation of the lattice vectors with 
respect to the Cartesian reference. See details.}

\item{x}{a 3x3 numeric matrix containing by column the coordinates of the 
lattice vectors.}
}
\value{
\code{cell2mat} returns a 3x3 numeric matrix containing by column the
  coordinates of the lattice vectors.
  
  \code{mat2cell} returns a data.frame containing by column the lattice
  parameters a, b, c, alpha, beta, gamma.
}
\description{
Convert lattice parameters into lattice vectors and vice versa.
}
\note{
Lengths are in Angstrom and angles in degree.
}

