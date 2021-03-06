% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/areColors.R
\name{areColors}
\alias{areColors}
\title{Are the elements of a vector valid R colors?}
\usage{
areColors(x)
}
\arguments{
\item{x}{an atomic vector to be tested.}
}
\value{
a logical vector indicating if the elements of \code{x} are valid
  colors.
}
\description{
Check if the elements of a vector are valid R colors.
}
\examples{
areColors(c(1, "1", NA, "black", "#000000"))

}
\seealso{
\code{\link[grDevices]{colors}}
}
\keyword{utilities}

