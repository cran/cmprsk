\name{crr}
\alias{crr}
\title{
Competing Risks Regression
}
\description{
regression modeling of subdistribution functions in competing risks
}
\usage{
crr(ftime, fstatus, cov1, cov2, tf, cengroup, failcode=1, cencode=0,
subset, na.action=na.omit, gtol=1e-06, maxiter=10, init)
}
\arguments{
\item{ftime}{
vector of failure/censoring times
}
\item{fstatus}{
vector with a unique code for each failure type and a separate code for
censored observations 
}
\item{cov1}{
matrix (nobs x ncovs) of fixed covariates (either cov1, cov2, or both
are required)
}
\item{cov2}{
matrix of covariates that will be multiplied by functions of time; 
if used, often these covariates would also appear in cov1
to give a prop hazards effect plus a time interaction
}
\item{tf}{
functions of time.  A function that takes a vector of times as
an argument and returns a matrix whose jth column is the value of 
the time function corresponding to the jth column of cov2 evaluated
at the input time vector.  At time \code{tk}, the
model includes the term \code{cov2[,j]*tf(tk)[,j]} as a covariate.
}
\item{cengroup}{
vector with different values for each group with 
a distinct censoring distribution (the censoring distribution
is estimated separately within these groups).  All data in one group, if
missing. 
}
\item{failcode}{
code of fstatus that denotes the failure type of interest
}
\item{cencode}{
code of fstatus that denotes censored observations
}
\item{subset}{
  a logical vector specifying a subset of cases to include in the
  analysis
}
\item{na.action}{
  a function specifying the action to take for any cases missing any of
  ftime, fstatus, cov1, cov2, cengroup, or subset.
}
\item{gtol}{
  iteration stops when a function of the gradient is \code{< gtol}
}
\item{maxiter}{
maximum number of iterations in Newton algorithm (0 computes
scores and var at \code{init}, but performs no iterations)
}
\item{init}{
  initial values of regression parameters (default=all 0)
}
}
\value{
  Returns a list of class crr, with components
  \item{$coef}{the estimated regression coefficients}
  \item{$loglik}{log pseudo-liklihood evaluated at coef}
\item{$score}{derivitives of the log pseudo-likelihood evaluated at coef}
\item{$inf}{-second derivatives of the log pseudo-likelihood}
\item{$var}{estimated variance covariance matrix of coef}
\item{$res}{matrix of residuals giving the
contribution to each score (columns) at each unique failure time
(rows)}
\item{$uftime}{vector of unique failure times}
\item{$bfitj}{jumps in the Breslow-type estimate of the underlying
  sub-distribution cumulative hazard (used by predict.crr())}
\item{$tfs}{the tfs matrix (output of tf(), if used)}
\item{$converged}{TRUE if the iterative algorithm converged.}
}
\details{
Fits the 'proportional subdistribution hazards' regression model
described in Fine and Gray (1999).  This model directly assesses the
effect of covariates on the subdistribution of a particular type of
failure in a competing risks setting.  The method implemented here is
described in the paper as the weighted estimating equation.

The basic model assumes the subdistribution with covariates z is a
constant shift on the complementary log log scale from a baseline
subdistribution function.  This can be generalized by including
interactions of z with functions of time to allow the magnitude of the
shift to change with follow-up time, through the cov2 and tfs
arguments.  For example, if z is a vector of covariate values, and uft
is a vector containing the unique failure times for failures of the
type of interest (sorted in ascending order), then the coefficients a,
b and c in the quadratic (in time) model
\eqn{az+bzt+zt^2}{a*z+b*z*t+c*z*t*t} can be fit
by specifying \code{cov1=z}, \code{cov2=cbind(z,z)},
\code{tf=function(uft) cbind(uft,uft*uft)}. 

This function uses an estimate of the survivor function of the
censoring distribution to reweight contributions to the risk sets for
failures from competing causes.  In a generalization of the methodology
in the paper, the censoring distribution can be estimated separately
within strata defined by the cengroup argument.  If the censoring
distribution is different within groups defined by covariates in the
model, then validity of the method requires using separate estimates of
the censoring distribution within those groups.

The residuals returned are analogous to the Schoenfeld residuals in
ordinary survival models.  Plotting the jth column of res against the
vector of unique failure times checks for lack of fit over time in the
corresponding covariate (column of cov1).
}
\references{
Fine JP and Gray RJ (1999) A proportional hazards model for the
subdistribution of a competing risk.  JASA 94:496-509.
}
\seealso{
\code{\link{predict.crr}} \code{\link{print.crr}} \code{\link{plot.predict.crr}}
}
\examples{
# simulated data to test 
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(600),nrow=200)
print(z <- crr(ftime,fstatus,cov))
# convergence:  TRUE 
# coefficients:
# [1]  0.8470 -0.3473 -0.6913
# standard errors:
# [1] 0.4027 0.4402 0.4916
# two-sided p-values:
# [1] 0.035 0.430 0.160
# for failures of type 2
z.p <- predict(z,rbind(c(.1,.5,.8),c(.1,.5,.2)))
plot(z.p,lty=1,color=2:3)
crr(ftime,fstatus,cov,failcode=2)
# convergence:  TRUE 
# coefficients:
# [1] -0.4638  0.4896  0.2331
# standard errors:
# [1] 0.3896 0.4540 0.4256
# two-sided p-values:
# [1] 0.23 0.28 0.58
# quadratic in time for first cov
crr(ftime,fstatus,cov,cbind(cov[,1],cov[,1]),function(Uft) cbind(Uft,Uft^2))
# convergence:  TRUE 
# coefficients:
# [1]  1.3360 -0.3484 -0.6852 -1.1210  0.3897
# standard errors:
# [1] 0.8657 0.4405 0.4909 1.6860 0.6271
# two-sided p-values:
# [1] 0.12 0.43 0.16 0.51 0.53
#last 2 values are the coefs of the cov[,1]*t and cov[,1]*t^2 terms
#additional examples in test.R
}
\keyword{survival}
