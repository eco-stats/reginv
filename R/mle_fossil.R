#' Confidence Interval for Extinction Time using maximum likelihood (asymptotic)
#'
#' Estimates a confidence interval for extinction time using maximum likelihood techniques -- either a Wald interval is used or
#' a profile likelihood technique is used, inverting the likelihood ratio statistic comparing to the appropriate quantile for the chi-squared distribution.
#' These methods do not have good small sample properties, especially when the measurement error is small, \code{\link{reginv_fossil}} is preferred.
#'
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param alpha Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals).
#'  If \code{alpha=NULL}, returns a maximum likelihood estimator only.
#' @param wald logical; FALSE (default) uses profile likelihood, comparing the likelihood ratio statistic to a quantile from the chi-squared distribution,
#'  TRUE uses a Wald interval (which has very poor performance for small sample sizes!)
#' @param ... logical; FALSE (default) uses profile likelihood, comparing the likelihood ratio statistic to a quantile from the chi-squared distribution,
#'  TRUE uses a Wald interval (which has very poor performance for small sample sizes!)
#'
#' @details 
#' Given a vector of fossil ages \code{ages} and corresponding measurement error standard deviations \code{sd}, and an upper limit \code{K} for the possible age of a fossil
#' that could be included in this dataset, our procedure then assumes that:
#' \itemize{
#' \item For a fossil of a given age, measurement error estimating its age is normally distributed with mean zero and standard deviation as provided
#' \item That fossil dates are uniformly distributed over the interval of allowable dates
#' \item All observed fossil dates are less than \code{K}
#' }
#' We then estimate extinction time using the maximum likelihood estimator, and construct a confidence interval either using profile likelihood
#' or a Wald interval. The profile likelihood approach finds the set of all extinction times \code{theta} such that a likelihood ratio test that the true
#' extinction time were \code{theta} would not be significant at level \code{alpha}, when comparing the likelihood ratio to a chi-square distribution with one degree of freedom.
#' The Wald interval assumes the maximum likelihood estimator is normally distributed and constructs a symmetric interval using the observed information matrix.
#' The Wald interval performs poorly unless measurement error and sample size are both large and is not recommended (except perhaps to obtain starting estimates for other algorithms).
#' 
#' It is assumed that \code{ages} has been specified with smaller values representing more recent specimens, for example, \code{ages} could be specified in years before present.
#' If there is interest in estimating speciation or invasion time, data would only need to be reordered so that smaller values represent older specimens. 
#' 
#' @return This function returns an object of class "reginv" with the following components:
#'
#'  \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}. (If \code{q} was not specified as input, this defaults to the lower limit for a \code{100(1-alpha)}\% confidence interval, a point estimate at \code{q=0.5} ("best estimate" of extinction time), and an upper limit for a \code{100(1-alpha)}\% confidence interval.)}
#'  \item{q}{ the vector of quantiles used in estimation.}
#'  \item{call }{ the function call}
#' @export
#' @examples
#' ages = rfossil(20, 10000, K=25000, sd=1000) #simulating some random data
#' # for the maximum likelihood estimate and a 95% CI
#' mle_fossil(ages, sd=1000, K=25000, alpha=0.05) 
#' # get the MLE only
#' mle_fossil(ages, sd=1000, K=25000, alpha=NULL) 

mle_fossil = function(ages, sd, K, alpha=0.05, wald=FALSE, ...)
{  
  thetaMLE = getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K)
  SE = 1/sqrt(-thetaMLE$hessian)
  if(is.null(alpha))
  {
    result = list(theta=thetaMLE$par, se=SE)
  }
  else
  {
    if(wald==TRUE)
    {
      lo = list( root=thetaMLE$par - qnorm(1-alpha/2) * SE )
      hi = list( root=thetaMLE$par + qnorm(1-alpha/2) * SE )
    }
    else
    {
      lo = try( uniroot(fossil_LRT,thetaMLE$par*c(0.25,1),thetaMLE,alpha=alpha, ages=ages,sd=sd,K=K,extendInt="downX", ...) )
      if(inherits(lo,"try-error")) lo=list(root=thetaMLE$par)
      hi = try( uniroot(fossil_LRT,thetaMLE$par*c(1,1.25),thetaMLE,alpha=alpha, ages=ages,sd=sd,K=K,extendInt="upX", ...) )
      if(inherits(hi,"try-error")) hi=list(root=thetaMLE$par)
      result = list( theta=c(lower=lo$root,mle=thetaMLE$par,upper=hi$root), se=SE)
    }
    result$call <- match.call()
    class(result)="reginv"
  }
  return( result )
}

getThetaMLE = function(ages, theta=min(ages), sd, K )
{
  if(all(sd==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,fossil_LogLik,ages=ages,sd=sd,K=K,method="Brent",lower=min(theta*c(0,2)),upper=max(theta*c(0,2)),control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}

fossil_LogLik = function(theta,ages,sd,K)
{
  return(ll=sum(dfossil(x=ages,theta=theta,K=K,sd=sd,log.p=TRUE)))
}

fossil_LRT = function(theta0,thetaMLE,ages, sd, K, alpha=0.05)
{
  ll0=fossil_LogLik(theta0,ages,sd,K)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
}
