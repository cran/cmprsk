\name{predict.crr}
\alias{predict.crr}
\title{
Estimate subdistribution functions from crr output
}
\description{
predict method for crr
}
\usage{
predict.crr(z, cov1, cov2)
}
\arguments{
\item{z}{
output from crr
}
\item{cov1, cov2}{
each row of cov1 and cov2 is a set of covariate values where the
subdistribution should be estimated.  The columns of cov1 and cov2 must
be in the same order as in the original call to crr.  Each must be
given if present in the original call to crr.
}}
\value{
Returns a matrix with the unique type 1 failure times in the first
column, and the other columns giving the estimated subdistribution
function corresponding to the covariate combinations in the rows of cov1
and cov2, at each failure time (the value that the estimate jumps to at
that failure time).
}
\details{
Computes \eqn{1-\exp(-B(t))}{1-exp(-B(t))}, where \eqn{B(t)}{B(t)} is
the estimated cumulative 
sub-distribution hazard obtained for the specified covariate values,
obtained from the Breslow-type estimate of the underlying hazard and
the estimated regression coefficients.
}
\seealso{
\code{\link{crr}} \code{\link{plot.predict.crr}}
}
\keyword{survival}
% Converted by Sd2Rd version 1.10.
