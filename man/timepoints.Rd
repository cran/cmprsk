\name{timepoints}
\alias{timepoints}
\title{
Calculate Estimates at Specific Timepoints
}
\description{
Find values at specified timepoints from curves specified as all
corners of step functions.
}
\usage{
timepoints(w, times)
}
\arguments{
\item{w}{
a list containing the estimates, with points for all corners of the
step function.  (Usually created by cuminc.)
Each component in the list contains the estimate for a different group.
Each component has components giving times, function estimates, and
variances (see cuminc)
}
\item{times}{
vector of times where estimates are needed
}}
\value{
  A list with components
  \item{$est}{a matrix of estimates of the subdistributions with a row for
each component in \code{w} and a column for each time}
  \item{$var}{a matrix giving the corresponding variances.}
}
\seealso{
\code{\link{cuminc}}
}
\keyword{survival}
% Converted by Sd2Rd version 1.10.
