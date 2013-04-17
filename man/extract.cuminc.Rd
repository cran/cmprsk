\name{[.cuminc}
\alias{[.cuminc}
\title{
Subset method for lists of class cuminc
}
\description{
A subset method that preserves the class of objects of class cuminc,
allowing a subset of the curves to be selected.
}
\usage{
\method{[}{cuminc}(x,i,\ldots)
}
\arguments{
  \item{x}{object of class cuminc}
  \item{i}{elements to extract}
  \item{...}{not used}
}

\value{
A list with selected components of \code{x}, with the class set to
cuminc so cuminc methods can be applied.
}
\seealso{
\code{\link{cuminc}} \code{\link{plot.cuminc}} \code{\link{print.cuminc}}
}

\keyword{survival}
