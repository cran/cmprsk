\name{print.cuminc}
\alias{print.cuminc}
\title{Print cuminc objects}
\description{
A print method for objects of class cuminc (output from \code{cuminc()}).
}
\usage{
\method{print}{cuminc}(x, ntp=4, maxtime, \dots)
}
\arguments{
 \item{x}{an object of class cuminc}
 \item{ntp}{number of timepoints where estimates are printed}
 \item{maxtime}{the maximum timepoint where values are printed.  The
   default is the maximum time in the curves in \code{x}}
 \item{...}{additional arguments to \code{print()}}
}
\details{
Prints the test statistics and p-values (if present in \code{x}), and for each
estimated cumulative incidence curve prints its value and estimated
variance at a vector of times.  The times are chosen between 0 and
maxtime using the \code{pretty()} function.
}
\author{Robert Gray}

\seealso{ \code{\link{cuminc}} }

%\examples{
%#see help(cuminc)
%}
\keyword{survival }%-- one or more ...
