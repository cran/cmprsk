\name{cuminc}
\alias{cuminc}
\title{
Cumulative Incidence Analysis
}
\description{
Estimate cumulative incidence functions from competing risks
data and test equality across groups
}
\usage{
cuminc(ftime, fstatus, group, strata, rho=0, cencode=0,
subset, na.action=na.omit)
}
\arguments{
\item{ftime}{
failure time variable
}
\item{fstatus}{
variable with distinct codes for different causes of failure
and also a distinct code for censored observations
}
\item{group}{
estimates will calculated within groups given by distinct values of this
variable.  Tests will compare these groups.  If missing then treated as all
one group (no test statistics)
}
\item{strata}{
stratification variable.  Has no effect on estimates.  Tests will be
stratified on this variable.  (all data in 1 stratum, if missing)
}
\item{rho}{
Power of the weight function used in the tests.
}
\item{cencode}{
value of fstatus variable which indicates the failure time is censored.
}
\item{subset}{
  a logical vector specifying a subset of cases to include in the
  analysis
}
\item{na.action}{
  a function specifying the action to take for any cases missing any of
  ftime, fstatus, group, strata, or subset.
}
}
\value{
A list with components giving the subdistribution estimates for each
cause in each group, and a component \code{Tests} giving the test
statistics and p-values for comparing the subdistribution for each cause
across groups (if the 
number of groups is \eqn{>}{>}1).  The components giving the estimates
have names that are a combination 
of the group name and the cause code.  
These components are also lists, with components
\item{\code{time}}{ the times
where the estimates are calculated}
\item{\code{est}}{the estimated
sub-distribution functions.  These are step functions (all corners
of the steps given), so they can be plotted using ordinary lines() commands.
Estimates at particular times can be located using the timepoints()
function.}
\item{\code{var}}{the estimated variance of
  the estimates, which are estimates of the asymptotic
  variance of Aalen (1978).  }
}
\references{
Gray RJ (1988) A class of K-sample tests for comparing the cumulative
incidence of a competing risk, ANNALS OF STATISTICS, 16:1141-1154.


Kalbfleisch and Prentice (1980) THE ANALYSIS OF FAILURE TIME DATA, p 168-9.


Aalen, O. (1978) Nonparametric estimation of partial transition
probabilities in multiple decrement models, ANNALS OF STATISTICS,
6:534-545.
}
\author{Robert Gray}
\seealso{
\code{\link{plot.cuminc}} \code{\link{timepoints}} \code{\link{print.cuminc}}
}
\examples{
set.seed(2)
ss <- rexp(100)
gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('a','b','c'))
cc <- sample(0:2,100,replace=TRUE)
strt <- sample(1:2,100,replace=TRUE)
print(xx <- cuminc(ss,cc,gg,strt))
plot(xx,lty=1,color=1:6)
# see also test.R, test.out
}
\keyword{survival}
