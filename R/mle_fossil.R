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
#' @param df Numeric; degrees of freedom for the t-distribution used to model measurement error. Must be at least 2. Default (NULL) uses a Gaussian distribution.
#' @param alpha Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals).
#'  If \code{alpha=NULL}, returns a maximum likelihood estimator only.
#' @param q Numeric vector of values between 0 and 1, specifying the quantiles at which we want to solve for extinction time. Defaults to \code{c(alpha/2,1-alpha/2)},
#' which gives the limits of a 100(1-\code{alpha})\% confidence interval. If \code{q} is specified it overrides any input for \code{alpha}.
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
#' @return This function returns an object of class "mle_fossil" with the following components:
#'
#'  \item{theta}{ a maximum likelihood estimator of \code{theta}.}
#'  \item{se}{ the estimated standard error of the MLE.}
#'  \item{ci}{ a vector of confidence limits for \code{theta} at the chosen confidence levels.}
#'  \item{q}{ the vector of quantiles used in estimation (if applicable).}
#'  \item{call }{ the function call}
#' @export
#' @examples
#' ages = rfossil(20, 10000, K=25000, sd=1000) #simulating some random data
#' # for the maximum likelihood estimate and a 95% CI
#' mle_fossil(ages, sd=1000, K=25000, alpha=0.05) 
#' # get the MLE only
#' mle_fossil(ages, sd=1000, K=25000, alpha=NULL) 

mle_fossil = function(ages, sd, K, df=NULL, alpha=0.05, q=c(alpha/2,1-alpha/2), wald=FALSE, ...)
{  
  nSD = length(sd)
  n = length(ages)
  if(nSD==1) sd = rep(sd,n)
  if(nSD!=n & nSD!=1)
  {
    sd = rep(sd,length=n)
    warning("lengths of 'ages' and sd' are not equal so will extend 'sd' as needed.")
  }
  
  if(all(ages>K)) stop("'ages' need to be no larger than 'K'")
  if(any(ages>K))
  {
    warning("Some ages exceed 'K', these will be ignored")
    sd   = sd[ages<=K]
    ages = ages[ages<=K]
  }
  thetaMLE = getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=df)
  SE = 1/sqrt(-thetaMLE$hessian)
  
  if(is.null(alpha))
  {
    result = list(theta=thetaMLE$par, se=SE, call=match.call())
  }
  else
  {
    # set up result list.
    ci=rep(NA,length(q))
    if(is.null(names(q)))
      names(ci)=paste0("q=",q)
    else
      names(ci)=names(q)

    if(wald==TRUE)
    {
      ci = thetaMLE$par + qnorm(q) * SE
    }
    else
    {
      nQ=length(q)
      # set search limits so that we look above MLE if q>0.5 and below otherwise 
      q2Tail = 2*pmin(q,1-q)
      qLo = qHi = rep(1,nQ)
      qLo[q<=0.5] = 0.25
      qHi[q>=0.5] = 1.25
      # note LRT function is increasing for q>0.5 
      dir = rep("downX",nQ)
      dir[q>=0.5]="upX"
      for(iQ in 1:nQ)
      {
        thLim = try( uniroot(fossil_LRT,thetaMLE$par*c(qLo[iQ],qHi[iQ]),thetaMLE,alpha=q2Tail[iQ], ages=ages,sd=sd,K=K,df=df,extendInt=dir[iQ], ...) )
        if(inherits(thLim,"try-error"))
          ci[iQ] = thetaMLE$par
        else
          ci[iQ] = thLim$root
      }
    }
    result = list( theta=thetaMLE$par, ci=ci, se=SE, q=q, call=match.call())
  }
  class(result)="mle_fossil"
  return( result )
}

getThetaMLE = function(ages, theta=min(ages), sd, K, df )
{
  if(all(sd==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,fossil_LogLik,ages=ages,sd=sd,K=K,df=df,method="Brent",lower=-1/sqrt(.Machine$double.eps),upper=K,control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}

fossil_LogLik = function(theta,ages,sd,K,df)
{
  return(ll=sum(dfossil(x=ages,theta=theta,K=K,sd=sd,df=df,log=TRUE)))
}

fossil_LRT = function(theta0,thetaMLE,ages, sd, K, df, alpha=0.05)
{
  ll0=fossil_LogLik(theta0,ages,sd,K,df)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
}
