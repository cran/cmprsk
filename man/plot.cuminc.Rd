\name{plot.cuminc}
\alias{plot.cuminc}
\title{
Create Labeled Cumulative Incidence Plots
}
\description{
Plot method for cuminc.  Creates labeled line plots from appropriate
list input, for example, the output from \code{cuminc()}.
}
\usage{
plot.cuminc(x, main=" ", curvlab, ylim=c(0, 1), xlim, wh=2, xlab="Years", 
ylab="Probability", lty=1:length(x), color=1, \dots)
}
\arguments{
\item{x}{
a list, with each component representing one curve in the plot.  Each
component of \code{x} is itself a list whose first component gives the x values
and 2nd component the y values to be plotted.  Although written for
cumulative incidence curves, can in principle be used for any set of lines.
}
\item{main}{
the main title for the plot.
}
\item{curvlab}{
Curve labels for the plot.  Default is \code{names(x)}, or if that is missing,
\code{1:nc}, where \code{nc} is the number of curves in \code{x}.
}
\item{ylim}{
yaxis limits for plot
}
\item{xlim}{
xaxis limits for plot (default is 0 to the largest time in any of the
curves)
}
\item{wh}{
  if a vector of length 2, then the upper right coordinates of the
  legend; otherwise the legend is placed in the upper right corner of
  the plot
}
\item{xlab}{
X axis label 
}
\item{ylab}{
 y axis label
}
\item{lty}{
vector of line types.  Default \code{1:nc} (\code{nc} is the number of
curves in \code{x}).  For color displays, \code{lty=1, color=1:nc}, might
be more appropriate.  If \code{length(lty)<nc}, then \code{lty[1]} is
used for all. 
}
\item{color}{
  vector of colors.  If \code{length(color)<nc}, then the \code{color[1]} is
  used for all.
  }
\item{...}{
additional arguments passed to the initial call of the plot function.
}}
\value{
No value is returned.  
}
\seealso{ \code{\link{cuminc}} }
%\examples{
%#see help(cuminc)
%}
\keyword{survival}
\keyword{hplot}
% Converted by Sd2Rd version 1.10.
