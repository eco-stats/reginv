#' Confidence Interval for Extinction Time using maximum likelihood (asymptotic)
#'
#' Estimates a confidence interval for extinction time using maximum likelihood techniques -- either a Wald interval is used or
#' a profile likelihood technique is used, inverting the likelihood ratio statistic comparing to the appropriate quantile for the chi-squared distribution.
#' These methods do not have good small sample properties, especially when the measurement error is small, \code{\link{reginv_cutt}} is preferred in these settings.
#'
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param df Numeric; degrees of freedom for the t-distribution used to model measurement error. Set to NULL to estimate degrees of freedom from the data. If a number, must be greater than 1. Default (Inf) uses a Gaussian distribution.
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
#' @seealso rcutt
#' @return This function returns an object of class "mle_cutt" with the following components:
#'
#'  \item{theta}{ a maximum likelihood estimator of \code{theta}.}
#'  \item{se}{ the estimated standard error of the MLE.}
#'  \item{ci}{ a vector of confidence limits for \code{theta} at the chosen confidence levels.}
#'  \item{q}{ the vector of quantiles used in estimation (if applicable).}
#'  \item{call }{ the function call}
#' @export
#' @examples
#' ages = rcutt(20, 10000, K=25000, sd=1000) #simulating some random data
#' # for the maximum likelihood estimate and a 95% CI
#' mle_cutt(ages, sd=1000, K=25000, alpha=0.05) 
#' # get the MLE only
#' mle_cutt(ages, sd=1000, K=25000, alpha=NULL) 

mle_cutt = function(ages, sd, K, df=Inf, alpha=0.05, q=c(alpha/2,1-alpha/2), wald=FALSE, ...)
{  
  dfMin=2
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
  if(all(sd==0) & is.null(df)) df=Inf # if sds are zero then df is irrelevant
  if(is.null(df))
  {
    dfKnown = FALSE
    mles = getJointMLE(ages=ages, theta=min(ages), sd=sd, K=K, dfMin=dfMin)
    thetaMLE = list(par=mles$par[1])
    dfOut = 1/mles$par[2]
    vr = try( solve(-mles$hessian) )
    if(inherits(vr,"try-error")) vr=matrix(NaN,2,2)
    SE = if(is.nan(vr[1,1])) 0 else sqrt(vr[1,1])
  }
  else
  {
    dfKnown = TRUE
    dfOut=df
    thetaMLE = getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=df)
    SE = 1/sqrt(-thetaMLE$hessian)
  }

  if(is.null(alpha))
  {
    result = list(theta=thetaMLE$par, se=SE, df=dfOut, call=match.call())
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
      nQ = length(q)
      q2Tail      = 2*pmin(q,1-q)
      # set search limits so that we look above MLE if q>0.5 and below otherwise 
      is_SE_bad   = is.nan(SE) | is.infinite(SE) | SE==0
      searchLim   = ifelse( is_SE_bad, IQR(ages)*0.5, SE*5 )
      qLo = qHi   = rep(thetaMLE$par,nQ)
      qLo[q<=0.5] = thetaMLE$par-searchLim
      qHi[q>=0.5] = min(thetaMLE$par+searchLim,K)
      # note LRT function is increasing for q>0.5 
      dir         = rep("downX",nQ)
      dir[q>=0.5] ="upX"
      for(iQ in 1:nQ)
      {
        if(dfKnown)
          thLim = try( uniroot(cutt_LRT, c(qLo[iQ],qHi[iQ]), thetaMLE, alpha=q2Tail[iQ], ages=ages, sd=sd, K=K, df=df, extendInt=dir[iQ], ...) )
        else
          thLim = try( uniroot(cutt_LRTprofile, c(qLo[iQ],qHi[iQ]), mles, alpha=q2Tail[iQ], ages=ages, sd=sd, K=K, dfMin=dfMin, extendInt=dir[iQ], ...) )
        if(inherits(thLim,"try-error"))
          ci[iQ] = thetaMLE$par
        else
          ci[iQ] = thLim$root
      }
    }
    result = list( theta=thetaMLE$par, ci=ci, se=SE, q=q, df=dfOut, call=match.call())
  }
  class(result)="mle_cutt"
  return( result )
}

getThetaMLE = function(ages, theta=min(ages), sd, K, df )
{
  if(all(sd==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,cutt_LogLik,ages=ages,sd=sd,K=K,df=df,method="Brent",lower=-1/sqrt(.Machine$double.eps),upper=K,control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}

getJointMLE = function(ages, theta=min(ages), sd, K, df=Inf, nIter=10, tol=1.e-5, dfMin=2 )
{
  if(all(sd==0))
    MLE = list( par=c(min(ages),Inf), value=length(ages)*log(1/(K-min(ages))), hessian=matrix(-Inf,2,2) )
  else
  {
    iIter = 1
    cond  = FALSE
    MLE   = optim(theta, cutt_LogLik, 
                 ages=ages, sd=sd, K=K, df=df, method="Brent",
                 lower=-1/sqrt(.Machine$double.eps), upper=K, control=list(trace=TRUE,fnscale=-1))
    while(cond==FALSE)
    {
      pre   = MLE
      df    = getDF( ages, theta=MLE$par, sd=sd, K=K, dfInvInit=1/df, dfMin=dfMin )$par
      MLE   = optim( pre$par, cutt_LogLik, 
                   ages=ages, sd=sd, K=K, df=df, method="Brent",
                   lower=-1/sqrt(.Machine$double.eps), upper=K, control=list(trace=TRUE,fnscale=-1) )
      eps   = abs(pre$value-MLE$value)
      cond  = eps<tol | iIter>=nIter
      iIter = iIter+1
    }
    if(eps>tol) MLE$convergence = 1 else MLE$convergence = 0
    MLE$par = c( MLE$par, 1/df )
    MLE$hessian = optimHess( MLE$par, cutt_LogLikJoint, ages=ages, sd=sd, K=K, dfMin=dfMin, control=list(trace=TRUE,fnscale=-1) )
  }
  return(MLE)
}

getDF = function( ages, theta, sd, K, dfInvInit=0, dfMin=2 )
{
  if(all(sd==0))
    res = list( par=Inf, value=length(ages)*log(1/(K-min(ages))) )
  else
  {
    dfInv = optim( dfInvInit, cutt_LogLikT, ages=ages, sd=sd, theta=theta, K=K, dfMin=dfMin, method="Brent", 
                   lower=0.005, upper=1/dfMin-sqrt(.Machine$double.eps), control=list(trace=TRUE,fnscale=-1, maxit=10) )
    res = list(par=1/dfInv$par,value=dfInv$value)
  }
  return(res)
}

cutt_LogLik = function(theta,ages,sd,K,df)
{
  return(ll=sum(dcutt(x=ages,theta=theta,K=K,sd=sd,df=df,log=TRUE)))
}

cutt_LRT = function(theta0,thetaMLE,ages, sd, K, df, alpha=0.05)
{
  ll0 = cutt_LogLik(theta0,ages,sd,K,df)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
}

cutt_LogLikJoint = function(params,ages,sd,K,dfMin=2) #parameters are (theta, dfInv) 
{
  if(params[2]>=1/dfMin)
    ll=sum(dcutt(x=ages,theta=params[1],K=K,sd=sd,df=dfMin+sqrt(.Machine$double.eps),log=TRUE))*(1+params[2]-1/dfMin) # game it away from df=2
  else
    ll=sum(dcutt(x=ages,theta=params[1],K=K,sd=sd,df=1/params[2],log=TRUE))
  return(ll)
}

cutt_LogLikT = function(dfInv,ages,sd,theta,K,dfMin=2)
{
  if(dfInv>=1/dfMin)
    ll=sum(dcutt(x=ages,theta=theta,K=K,sd=sd,df=2+sqrt(.Machine$double.eps),log=TRUE))*(1+dfInv-1/dfMin) # game it away from df=2
  else
    ll=sum(dcutt(x=ages,theta=theta,K=K,sd=sd,df=1/dfInv,log=TRUE))
  return(ll)
}

cutt_LRTprofile = function(theta0, mles, ages, sd, K, alpha=0.05, dfMin=2)
{
  ll0 = getDF(ages,theta=theta0, sd=sd, K=K, dfInvInit = mles$par[2], dfMin=dfMin)$value
  return(-2*(ll0-mles$value)-qchisq(1-alpha,1))
}