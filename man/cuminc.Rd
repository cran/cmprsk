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
cause in each group, and a component giving the test results (if the
number of groups is >1).  The components giving the estimates have names
that are a combination 
of the group name and the cause code.  The estimate for each group and
cause combination is also a list.
\item{$<groupcause>}{
\item{\code{$<groupcause>$time}}{ the times
where the estimates are calculated}
\item{\code{$<groupcause>$est}}{the estimated
sub-distribution functions.  These are step functions (all corners
of the steps given), so they can be plotted using ordinary lines() commands.
Estimates at particular times can be located using the timepoints()
function.}
\item{\code{$<groupcause>$var}}{the estimated variance of
  the estimates, which are estimates of the asymptotic
  variance of Aalen (1978).  }
}

\item{$Tests}{
  If the number of groups is \eqn{>}{>}1, there is also a component at
  the top level giving the test
statistics and p-values for comparing the subdistribution for each
cause across groups.  The test statistics are described in Gray (1988).
}
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
#Tests:
#      stat        pv
#1 3.393977 0.1832345
#2 1.989511 0.3698139
#Estimates and Variances:
#$est
#            1         2         3         4         5
#a 1 0.1311269 0.2699184 0.3420625 0.3420625 0.3420625
#b 1 0.2176471 0.2615686        NA        NA        NA
#c 1 0.3816280 0.4889137 0.5723581        NA        NA
#a 2 0.2257601 0.2972171 0.4415053 0.4415053 0.4415053
#b 2 0.3117647 0.5878431        NA        NA        NA
#c 2 0.2160508 0.2607532 0.2607532        NA        NA
#
#$var
#              1           2           3          4          5
#a 1 0.003922836 0.009113464 0.012507959 0.01250796 0.01250796
#b 1 0.005528913 0.006916070          NA         NA         NA
#c 1 0.009854126 0.012351948 0.015575726         NA         NA
#a 2 0.005947601 0.007356628 0.014479506 0.01447951 0.01447951
#b 2 0.007033454 0.010927780          NA         NA         NA
#c 2 0.006497566 0.007867543 0.007867543         NA         NA
plot(xx,lty=1,color=1:6)
# see also test.R, test.out
}
\keyword{survival}
