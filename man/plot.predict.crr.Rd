\name{plot.predict.crr}
\alias{plot.predict.crr}
\title{
Plot estimated subdistribution functions
}
\description{
plot method for \code{predict.crr}
}
\usage{
plot.predict.crr(z, lty=1:(ncol(z)-1), color=1,  ylim=c(0, max(z[, -1])),
xmin=0, xmax=max(z[, 1]), \dots)
}
\arguments{
\item{z}{
Output from \code{predict.crr}
}
\item{lty}{
vector of line types. If length is \eqn{<}{<} \# curves, then
\code{lty[1]} is used for all. 
}
\item{color}{
vector of line colors.  If length is \eqn{<}{<} \# curves, then
\code{color[1]} is used for all. 
}
\item{ylim}{
range of y-axis (vector of length two)
}
\item{xmin}{
lower limit of x-axis (often 0, the default)
}
\item{xmax}{
upper limit of x-axis
}
\item{...}{
Other arguments to plot
}}
\section{Side Effects}{
plots the subdistribution functions estimated by \code{predict.crr}, by
default using a different line type for each curve
}
\seealso{
\code{\link{crr}} \code{\link{predict.crr}}
}
\keyword{survival}
% Converted by Sd2Rd version 1.10.
