\name{print.crr}
\alias{print.crr}
\title{
prints summary of a crr object
}
\description{
print method for crr objects
}
\usage{
print.crr(x, \dots)
}
\arguments{
\item{x}{
crr object (output from \code{crr()})
}
\item{...}{additional arguments to \code{print()}}
}
\details{
prints the convergence status, the estimated coefficients, the
estimated standard errors, and the two-sided p-values for the test of
the individual coefficients equal to 0. (If convergence is false
everything else may be meaningless.)
}
\seealso{
\code{\link{crr}}
}
\keyword{survival}
% Converted by Sd2Rd version 1.10.
