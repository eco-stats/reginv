#' Confidence Interval for Extinction Time using regression inversion
#'
#' Estimates a confidence interval for extinction time using regression inversion, which finds a range of values for extinction time that
#' are plausible considering the value of the maximum likelihood estimator. This is done assuming fossil dates come from a
#' a compound uniform-truncated T distribution, which accounts for sampling error (the fact that the most recent fossil date is not necessarily
#' the most recent time that the species was extant) and measurement error (error dating fossils).
#' Usually takes a while to run (but less than a minute), because it is doing a lot of computation behind the scenes. 
#'
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param df Numeric; degrees of freedom for the t-distribution used to model measurement error. Must be greater than 2. Uses a Gaussian distribution if \code{df=Inf}.
#' Default (\code{NULL}) estimates \code{df} from the data.
#' @param alpha Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals)
#' @param q Numeric vector of values between 0 and 1, specifying the quantiles at which we want to solve for extinction time. Defaults to \code{c(alpha/2,0.5,1-alpha/2)},
#' which gives the limits of a 100(1-\code{alpha})\% confidence interval and a point estimate obtained by solving at 0.5. If \code{q} is specified it overrides any input for \code{alpha}.
#' @param paramInits A numeric vector of initial values for extinction time to use in simulation. If \code{NULL} then these will be 100 values evenly distributed within 5 SEs of the estimate from \code{\link{mle_cutt}}.
#' @param iterMax Maximum number of simulated datasets to use to estimate extinction time (default 500).
#' @param method Regression method to use in estimating how MLE quantiles vary with extinction time, using a dataset of simulated MLEs and the extinction times at which they were simulated.
#'  \code{method='rq'} (default) uses linear quantile regression, \code{method='rq2'} uses quadratic quantile regression, \code{method='wrq'} uses linear quantile regression but down-weighting 
#'  high influence points, \code{method='prob'} uses linear probit regression on an indicator variable for whether or not simulated MLEs exceeds the observed MLE.
#'  \code{method='lm'} uses a linear regression, which is appropriate for bias correction rather than confidence interval estimation. In this case \code{q} and \code{alpha} are ignored and only a point estimate is returned.
#'
#' @details 
#' Given a vector of fossil ages \code{ages} and corresponding measurement error standard deviations \code{sd}, and an upper limit \code{K} for the possible age of a fossil
#' that could be included in this dataset, our procedure then assumes:
#' \itemize{
#' \item That, for a fossil of a given age, measurement error estimating its age follows a distribution that is Student's T multiplied by the provided standard deviation, truncated such that observed age is less than \code{K}
#' \item That, given a value for measurement error, fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until the value such that observed age is \code{K})
#' }
#' We then estimate extinction time by inversion of the maximum likelihood estimator, that is, we find the estimate of extinction time \eqn{\theta}{t}
#' such that if this were the true extinction time then the chance of seeing a MLE less than the one obtained from sample \code{ages} is equal to \code{q}. The probability
#' is estimated by regression, hence we call our procedure regression inversion. If \code{alpha} is specified this function will return three values:
#' the lower limit of the \code{100*(1-alpha)}\% confidence interval (solving at \code{q=alpha/2}), a point estimator for extinction time (solving at \code{q=0.5}), and an upper limit for the
#' confidence interval (solving at \code{q=1-alpha/2}). If a vector \code{q} is specified as input then the function solves for this vector instead. 
#' 
#' Note that because we are using simulation to estimate parameters and their confidence intervals, we will get slightly different answers on each run. 
#' 
#' It is assumed that \code{ages} has been specified with smaller values representing more recent specimens, for example, \code{ages} could be specified in years before present.
#' If there is interest in estimating speciation or invasion time, data would only need to be reordered so that smaller values represent older specimens. 
#' 
#' @return This function returns an object of class \code{est_cutt} with the following components:
#'
#'  \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}. (If \code{q} was not specified as input, this defaults to a bias-corrected maximum likelihood estimate \code{point}, and together with confidence interval limits \code{lo} and \code{hi}.)}
#'  \item{q}{ the vector of quantiles used in estimation.}
#'  \item{error}{ Monte Carlo standard error estimating each of these values, as estimated from the regression. }
#'  \item{se}{ the estimated standard error of the MLE.}
#'  \item{iter}{ the number of iterations taken to estimate this value.}
#'  \item{converged}{ whether or not this converged, at the specified tolerance.}
#'  \item{df}{ the estimated degrees of freedom of the Student's T distribution for measurement error.}
#'  \item{method}{ \code{reginv}.}
#'  \item{data}{ a data frame with the data used in analysis, in columns labelled as \code{ages} and \code{sd}.}
#'  \item{call }{ the function call.}
#' @seealso est_cutt, mle_cutt, cutt
#' @export
#' @examples
#' ages = rcutt(20, 10000, K=25000, sd=500) #simulating some random data
#' 
#' # for a point estimate together with a 95% CI (only 200 iterations used so it runs quickly)
#' reginv_cutt(ages=ages, sd=500, K=25000, alpha=0.05, iterMax=200) 
#' 
#' # compare to estimates using asymptotic likelihood inference, which tend to
#' # be narrower and have poorer coverage (they miss the true value too often
#' # when n or sd is small):
#' mle_cutt(ages=ages, sd=500, K=25000, alpha=0.05) 

reginv_cutt = function(ages, sd, K, df=NULL, alpha=0.05, q=c(lo=alpha/2,point=0.5,hi=1-alpha/2), paramInits=NULL, iterMax=1000, method="rq")
{  
  
  # if method="lm" we just want a point estimate
  if(method=="lm")
  {
    q=0.5
    names(q)="theta"
  }

  # set up result list.
  theta = rep(NA,length(q))
  if(is.null(names(q)))
    names(theta)=paste0("q=",q)
  else
    names(theta)=names(q)
  result=list(theta=theta,q=q,error=theta,iter=theta,converged=theta)
  
  # if paramInits provided, get data and test stats for initial parameters just once
  n <- length(ages)
  stats=data.frame(theta=NULL,T=NULL,thetaEst=NULL)
  if(is.null(paramInits)==FALSE)
  {
    for(i in 1:length(paramInits) )
    {
      newDat = simFn_cutt(theta=paramInits[i],K=K,sd=sd,df=df,n=n,dat=ages)
      stats = rbind(stats, c(theta=paramInits[i], T = getThMLE(newDat, K=K,sd=sd,df=df,n=n), thetaEst=paramInits[i]) )
    }      
    names(stats)=c("theta","T","thetaEst")
  }

  # get results for each value of q
  for (iQ in 1:length(q)) {

    # get paramInits, if not provided
    if(is.null(paramInits))
    {
      ft.mle = mle_cutt(ages=ages, sd=sd, K=K, df=df, q=q[iQ])
      result$se = ft.mle$se
      result$df = ft.mle$df
      is_SE_bad = is.nan(ft.mle$se) | is.infinite(ft.mle$se) | ft.mle$se==0
      stepSize = ifelse( is_SE_bad, IQR(ages)*0.1, ft.mle$se )
      if(ft.mle$theta[1]+5*stepSize>K) #make sure paramInits don't exceed K
        paramInits = seq( ft.mle$theta[1] -5* stepSize, K, length=100) 
      else
        paramInits = ft.mle$theta[1] + stepSize*seq(-5,5,length=100) 
      stats = NULL
    }
    # call reginv
    resulti = reginv(data=ages,getT=getThMLE,simulateData=simFn_cutt,paramInits=paramInits,
                        q=q[iQ],iterMax=iterMax,K=K,sd=sd, df=df, n=n, dat=ages, method=method,stats=stats,eps=1.e-15)
    result$theta[[iQ]] = resulti$theta
    result$error[[iQ]] = resulti$error
    result$iter[[iQ]]  = resulti$iter
    result$converged[[iQ]]  = resulti$converged
  }
  result$data = data.frame(ages=ages, sd=sd)
  result$method="reginv"
  result$call <- match.call()
  class(result)="reginv"
  class(result)="est_cutt"
  return(result)
}


simFn_cutt = function (theta, K, sd, df, n=length(sd), dat)
{
  if(theta>K+mean(sd)) #trying to game it to push estimates away from K 
    W = rep(theta,n)
  else
  {
    if( is.null(df) )  df = getDF(dat, theta=theta, sd=sd, K=K)$par
    W = rcutt(theta, K, sd, df=df, n=n)
  }
  return(W)
}

getThMLE = function (ages, K, sd, df, n=length(ages), dat=NULL )
{
  if(all(ages>K+mean(sd)))
    theta = min(ages)
  else
    theta = mle_cutt(ages, sd, K, df=df, alpha=NULL )$theta
  return(theta)
}