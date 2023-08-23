#' Confidence Interval for Extinction Time using regression inversion
#'
#' Estimates a confidence interval for extinction time using regression inversion, which finds a range of values for extinction time that
#' are plausible considering the maximum likelihood estimator of the most recently observed fossil, accounting for sampling error (the fact that
#' the most recent fossil date is not necessarily the most recent time that the species was extant) and measurement error (error dating fossils).
#' Usually takes a few seconds to run, because it is doing a lot of computation behind the scenes. 
#'
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param df Numeric; degrees of freedom for the t-distribution used to model measurement error. Must be at least 2. Default (NULL) uses a Gaussian distribution.
#' @param alpha Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals)
#' @param q Numeric vector of values between 0 and 1, specifying the quantiles at which we want to solve for extinction time. Defaults to \code{c(alpha/2,0.5,1-alpha/2)},
#' which gives the limits of a 100(1-\code{alpha})\% confidence interval and a point estimate obtained by solving at 0.5. If \code{q} is specified it overrides any input for \code{alpha}.
#' @param paramInits A numeric vector of initial values for extinction time to use in simulation. If \code{NULL} then these will be 20 values evenly distributed within 5 SEs of the MLE.
#' @param iterMax Maximum number of simulated datasets to use to estimate extinction time (default 500).
#' @param method Regression method to use in estimating how MLE quantiles vary with extinction time, using a dataset of simulated MLEs and the extinction times at which they were simulated.
#'  \code{method='rq'} (default) uses linear quantile regression, \code{method='rq2'} uses quadratic quantile regression, \code{method='wrq'} uses linear quantile regression but down-weighting 
#'  high influence points, \code{method='prob'} uses linear probit regression on an indicator variable for whether or not simulated MLEs exceeds the observed MLE.
#'
#' @details 
#' Given a vector of fossil ages \code{ages} and corresponding measurement error standard deviations \code{sd}, and an upper limit \code{K} for the possible age of a fossil
#' that could be included in this dataset, our procedure then assumes:
#' \itemize{
#' \item That, for a fossil of a given age, measurement error estimating its age is normally distributed with mean zero and standard deviation as provided, truncated such that observed age is less than \code{K}
#' \item That, given a value for measurement error, fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until the value such that observed age is \code{K})
#' }
#' We then estimate extinction time by inversion of the maximum likelihood estimator, that is, we find the estimate of extinction time \eqn{\theta}{t}
#' such that if this were the true extinction time then the chance of seeing a MLE less than the one obtained from sample \code{ages} is equal to \code{q}. The probability
#' is estimated by regression, hence we call our procedure regression inversion. If \code{alpha} is specified this function will return three values:
#' the lower limit of the \code{100*(1-alpha)}\% confidence interval (solving at \code{q=alpha/2}), a point estimator for extinction time (solving at \code{q=0.5}), and an upper limit for the
#' confidence interval (solving at \code{q=1-alpha/2}). If a vector \code{q} is specified as input then the function solves for this vector instead. 
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
#' 
#' # for a point estimate together with a 95% CI
#' reginv_fossil(ages=ages, sd=1000, K=25000, alpha=0.05) 
#' 
#' # compare to estimates using asymptotic likelihood inference, which tend to
#' # be narrower and have poorer coverage (they miss the true value too often
#' # when n or sd is small):
#' mle_fossil(ages=ages, sd=1000, K=25000, alpha=0.05) 

reginv_fossil = function(ages, sd, K, df=NULL, alpha=0.05, q=c(alpha/2,0.5,1-alpha/2), paramInits=NULL, iterMax=500, method="rq")
{  
  # get paramInits, if not provided
  if(is.null(paramInits))
  {
    ft.mle = mle_fossil(ages=ages, sd=sd, K=K, df=df, alpha=NULL)
    stepSize = ifelse( any(sd==0), IQR(ages)*0.1, ft.mle$se )
    if(method=="prob") #give prob extra starting values to reduce chance of separation
      paramInits = ft.mle$theta + stepSize*seq(-5,5,length=100) 
    else
      paramInits = ft.mle$theta + stepSize*seq(-5,5,length=20)
  }
  
  # set up result list.
  theta = rep(NA,length(q))
  if(is.null(names(q)))
    names(theta)=paste0("q=",q)
  else
    names(theta)=names(q)
  
  n <- length(ages)
  result=list(theta=theta,q=q)
  for (iQ in 1:length(q)) {
    result$theta[[iQ]] = reginv(ages,getT=getThMLE,simulateData=simFn_fossil,paramInits=paramInits,
                        q=q[iQ],iterMax=iterMax,K=K,sd=sd, df=df, n=n, method=method)$theta
  }
  result$call <- match.call()
  class(result)="reginv"
  return(result)
}

simFn_fossil = function (theta, K, sd, df, n=length(sd))
{
  if(theta>K) #trying to game it to push estimates away from K 
    W = rep(theta,n)
  else
    W = rfossil(theta, K, sd, df=df, n=n)
  return(W)
}

getThMLE = function (ages, K, sd, df, n=length(ages) )
{
  if(all(ages>K))
    theta = min(ages)
  else
    theta = mle_fossil(ages, sd, K, df=df, alpha=NULL )$theta
  return(theta)
}