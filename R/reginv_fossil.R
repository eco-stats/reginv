#' Confidence Interval for Extinction Time using regression inversion
#'
#' Estimates a confidence interval for extinction time using regression inversion, which finds a range of values for extinction time that
#' are plausible considering the maximum likleihood estimator of the most recently observed fossil, accounting for sampling error (the fact that
#' the most recent fossil date is not necessarily the most recent time that the species was extant) and measurement error (error dating fossils).
#'
reginv_fossil = function(ages, sd, K, alpha, q=c(alpha/2,0.5,1-alpha/2), thetaInits=NULL, iterMax=500, method="rq")
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param alpha Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals)
#' @param q Numeric vector of values between 0 and 1, specifying the quantiles at which we want to solve for extinction time. Defaults to \code{c(alpha/2,0.5,1-alpha/2)},
#' which gives the limits of a 100(1-\code{alpha})\% confidence interval and a point estimate obtained by solving at 0.5. If \code{q} is specified it overrides any input for \code{alpha}.
#' @param thetaInits A numeric vector of initial values for extinction time to use in simulation. If \code{NULL} then these will be 20 values evenly distributed within 5 SEs of the MLE.
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
#' ages = runif(20, 10000, 25000) #simulating some random data
#' sd = runif(20, 50, 100)
#'
#' # for a point estimate plus 95% CI
#' reginv_fossil(ages=ages, sd=sd, K=22000, alpha=0.05) 

reginv_fossil = function(ages, sd, K, alpha, q=c(alpha/2,0.5,1-alpha/2), thetaInits=NULL, iterMax=500, method="rq")
{  
  # get thetaInits, if not provided
  if(is.null(thetaInits))
  {
    ft.mle = getThetaMLE(ages=ages, theta=min(ages), eps.sigma=sd, K=K)
    stepSize = max(1/sqrt(-ft.mle$hessian), IQR(ages)*0.1, na.rm=TRUE)
    thetaInits = ft.mle$par + stepSize*seq(-5,5,length=20)
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
    result$theta[[iQ]] = reginv(ages,getT=getThMLE,simulateData=simFn_fossil,thetaInits=thetaInits,
                        q=q[iQ],iterMax=iterMax,K=K,eps.sigma=sd, method=method)$theta
  }
  result$call <- match.call()
  class(result)="reginv"
  return(result)
}

simFn_fossil = function (theta, K, eps.sigma, n=length(eps.sigma))
{
  # Simulate fossils assuming:
  # - Gaussian measurement error (truncated at K-theta)
  # - Uniform deposition from theta to K-eps
  if(theta>K) #trying to game it to push estimates away from K 
    W=rep(theta,n)
  else
    W = rUNmod(theta, K, eps.sigma, n=length(eps.sigma))
  return(W)
}


rUNmod = function(theta, K, eps.sigma, n=length(eps.sigma), tol=sqrt(.Machine$double.eps), nIter=50)
{
  #ensure sds is the right length
  nSD = length(eps.sigma)
  if(nSD!=n)
  {
    eps.sigma = rep(eps.sigma,length=n)
    warning("length of 'eps.sigma' is not equal to 'n' so will extend 'eps.sigma' as needed.")
  }
  if(nSD==1)  eps.sigma = rep(eps.sigma,length=n)
  
  # compute C and u
  F.K = pnorm(K - theta, mean = 0, sd = eps.sigma)
  f.K = dnorm(K - theta, mean = 0, sd = eps.sigma)
  C = (K-theta) * F.K + eps.sigma^2 * f.K
  u = runif(n)
  uOnC = u * C
  
  # now get cracking finding eps
  epsOld = qnorm(u,mean=0,sd=eps.sigma) # starting estimate
  eps = quantU = rep(0,n) 
  quantUMax = pnorm(K-theta,sd=eps.sigma) # quantU can't get any larger than this or eps will be larger than K-theta
  iter=0
  isDiff = eps.sigma!=0 #only do the below when eps.sigma is non-zero
  while(any(isDiff) & iter<nIter)
  {
    quantU[isDiff] = ( uOnC[isDiff] - eps.sigma[isDiff]^2*dnorm(epsOld[isDiff], mean = 0, sd = eps.sigma[isDiff]) ) / (K-theta)
    quantU[isDiff] = pmax(quantU[isDiff], sqrt(tol)) # correction for wild estimates sending quantU negative
    quantU[isDiff] = pmin(quantU[isDiff], quantUMax[isDiff]-sqrt(tol)) # correction for wild estimates sending eps over K-theta
    eps[isDiff]    = qnorm(quantU[isDiff], mean=0, sd=eps.sigma[isDiff])
    epsDiff        = abs(eps-epsOld)
    isDiff         = epsDiff>tol
    epsOld[isDiff] = eps[isDiff]
    iter           = iter+1
  }
  
  if(iter==nIter)
  {
    for(iObs in which(isDiff))
    {
      epsTry = try( uniroot( pUNmodEps, interval=c(qnorm(sqrt(tol),mean=0,sd=eps.sigma[iObs]),K-theta),u=u[iObs],theta=theta,tol=tol,n=1,eps.sigma=eps.sigma[iObs],K=K,extendInt="upX") )
      if(inherits(epsTry,"try-error")==FALSE)
      {
        eps[iObs] = epsTry$root
        isDiff[iObs] = FALSE
      }
    }
    if(any(isDiff)) warning(paste0("non-convergence for ",sum(isDiff)," observations"))
  }
  X = runif(n,min=theta,max=K-eps)
  W = X + eps
  return(W)
}

pUNmodEps = function(eps,u,theta,K,eps.sigma,n=length(eps.sigma),tol=sqrt(.Machine$double.eps))
  # function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # get CDF denominator
  F.K    = pnorm(K - theta, mean = 0, sd = eps.sigma)
  f.K    = dnorm(K - theta, mean = 0, sd = eps.sigma)
  C      = (K-theta) * F.K + eps.sigma^2 * f.K
  
  # get CDF-u
  F.eps  = pnorm(eps, mean = 0, sd = eps.sigma)
  f.eps  = dnorm(eps, mean = 0, sd = eps.sigma)
  cdfEps = ( (K-theta)*F.eps + eps.sigma^2*f.eps ) / C - u
  return(cdfEps)
}

getThMLE = function (ages, theta=NULL, K, eps.sigma )
{
  if(is.null(theta))  theta=min(ages)
  if(all(ages>K))
    theta = min(ages)
  else
    theta = getThetaMLE(ages, theta=theta, eps.sigma, K )$par
  return(theta)
}

getThetaMLE = function(ages, theta=min(ages), eps.sigma, K )
{
  if(all(eps.sigma==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,UNmodLogLik,ages=ages,sds=eps.sigma,K=K,method="Brent",lower=min(theta*c(0,2)),upper=max(theta*c(0,2)),control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}

UNmodLogLik = function(theta,ages,sds,K)
{
  # get F(K-theta), f(K-theta), f(w-theta)
  lF.eps.w = pnorm(ages - theta, mean = 0, sd = sds, log.p=TRUE)
  F.eps.K = pnorm(K - theta, mean = 0, sd = sds)
  f.eps.K = dnorm(K - theta, mean = 0, sd = sds)
  C = (K-theta) * F.eps.K + sds^2 * f.eps.K
  dlUN = lF.eps.w - log(C)
  return(ll=sum(dlUN))
}