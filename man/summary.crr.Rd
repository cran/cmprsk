\name{summary.crr}
\alias{summary.crr}
\alias{print.summary.crr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for crr}
\description{
Generate and print summaries of crr output
}
\usage{
\method{summary}{crr}(object, conf.int = 0.95, digits =
max(options()$digits - 5, 2), ...)

\method{print}{summary.crr}(x, digits = max(options()$digits - 4, 3), ...)
}

\arguments{
  \item{object}{ An object of class crr (output from the crr function) }
  \item{conf.int}{the level for a two-sided confidence interval on the
    coeficients. Default is 0.95.   }
  \item{digits}{ In \code{summary.crr}, \code{digits} determines the number of
  significant digits retained in the p-values.  In
  \code{print.summary.crr},
  \code{digits} sets the values of the \code{digits} option for printing
  the output.}
  \item{\dots}{ Included for compatibility with the generic functions.  Not
    currently used. }
  \item{x}{An object of class summary.crr (output from the summary
  method for crr)}
}
\details{
  The summary method calculates the standard errors, subdistribution
  hazard ratios z-scores, p-values, and confidence intervals on the
  hazard ratios.  The print method prints a fairly standard format
  tabular summary of the results.

  The pseudo likelihood ratio test in the printed output is based on the
  difference in the objective function at the global null and at the
  final estimates.  Since this objective function is not a true
  likelihood, this test statistic is not asymptotically chi-square.
}
\value{\code{summary.crr} returns a list of class summary.crr, which
  contains components
  \item{call }{The call to crr}
  \item{converged }{TRUE if the iterative algorithm converged}
  \item{n }{The number of observations used in fitting the model}
  \item{n.missing }{The number of observations removed by \code{crr}
  from the input data due to missing values}
  \item{loglik }{The value of the negative of the objective function
  (the pseudo log likelihood at convergence}
  \item{coef }{A matrix giving the estimated coefficients, hazard
  ratios, standard errors, z-scores, and p-values }
  \item{conf.int }{A matrix giving the estimated hazard ratios, inverse
  hazard ratios and lower and upper confidence limits on the hazard ratios}
  \item{logtest }{Twice the difference in log pseudo likelihood values}
}
\author{The summary and print.summary methods were provided by Luca Scrucca}

\seealso{ \code{\link{crr}}}
\examples{
## see examples in the crr help file
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
