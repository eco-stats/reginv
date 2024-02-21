#' Confidence Interval for Extinction Time using CUTT distribution
#'
#' Estimates a point estimate and a confidence interval for extinction time using either regression inversion
#' or likelihood-based inference, given a vector of fossil ages together with their measurement errors. By default 
#' the decision about which method to use is determined by the size of measurement error relative to the length of 
#' the time interval for which we have fossils (likelihood-based inference used for larger measurement error). 
#' Extinction time is estimated using a maximum likelihood estimator for the compound uniform-truncated T 
#' distribution, which accounts for sampling error (the fact that the most recent fossil date is not necessarily
#' the most recent time that the species was extant) and measurement error (error dating fossils).
#' Regression inversion may take up to a minute to run, because it is doing a lot of computation behind the scenes. 
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
#' @param method Method used to construct point estimate and confidence intervals. Default is to use regression inversion (\code{\link{reginv_cutt}}) if
#' average measurement error is less than 10% of the difference between \code{K} and the minimum observed fossil age, otherwise,
#' likelihood-based inference is used via \code{\link{mle_cutt}}.
#' @param paramInits (for \code{method="reginv"} only) A numeric vector of initial values for extinction time to use in simulation. If \code{NULL} then these will be 100 values evenly distributed within 5 SEs of the estimate from \code{\link{mle_cutt}}.
#' @param iterMax (for \code{method="reginv"} only) Maximum number of simulated datasets to use to estimate extinction time (default 500).
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
#' is either estimated by regression inversion (if \code{method="reginv"}, using the \code{\link{reginv_cutt}} function) or using a chi-squared distribution (if \code{method="mle"}, 
#' using the \code{\link{mle_cutt}} function). Regression inversion is computationally intensive but works well in a broad range of settings, whereas likelihood inversion is quick but
#' is only reliable when measurement error is large, hence is the default in such settings (when average measurement error standard deviation is more than a tenth of the range of the data).
#' Note that because \code{method="reginv"} uses simulation to estimate parameters and their confidence intervals, we will get slightly different answers on different runs. 
#' 
#' 
#' If \code{alpha} is specified this function will return three values:
#' the lower limit of the \code{100*(1-alpha)}\% confidence interval (solving at \code{q=alpha/2}), a point estimator for extinction time (solving at \code{q=0.5}), and an upper limit for the
#' confidence interval (solving at \code{q=1-alpha/2}). For \code{method="reginv"}, a bias correction is used when finding the point estimator, also obtained via regression inversion. 
#' If a vector \code{q} is specified as input then the function solves for this vector instead. 
#' 
#' It is assumed that \code{ages} has been specified with smaller values representing more recent specimens, for example, \code{ages} could be specified in years before present.
#' If there is interest in estimating speciation or invasion time, data would only need to be reordered so that smaller values represent older specimens. 
#' 
#' @return This function returns an object of class "reginv" with the following components:
#'
#'  \item{theta}{ a maximum likelihood estimator of \code{theta}.}
#'  \item{se}{ the estimated standard error of the MLE.}
#'  \item{q}{ the vector of quantiles used in estimation (if applicable).}
#'  \item{df}{ the estimated degrees of freedom of the Student's T distribution for measurement error.}
#'  \item{data}{ a data frame with the data used in analysis, in columns labelled as \code{ages} and \code{sd}}
#'  \item{call }{ the function call}
#'  \item{method}{ the method used for estimation}
#' If \code{\link{reginv_cutt}} is used for estimation is used for estimation, we also get the following components:
#'  \item{error}{ Monte Carlo standard error estimating each of these values, as estimated from the regression. }
#'  \item{iter}{ the number of iterations taken to estimate this value}
#'  \item{converged}{ whether or not this converged, at the specified tolerance}

#' @seealso reginv_cutt, mle_cutt, cutt
#' @export
#' @examples
#' ages = rcutt(20, 10000, K=25000, sd=1000) #simulating some random data
#' 
#' # for a point estimate together with a 95% CI (only 200 iterations used so it runs quickly)
#' est_cutt(ages=ages, sd=500, K=25000, alpha=0.05, iterMax=200) 
#' 
#' # compare to estimates using asymptotic likelihood inference, which tend to
#' # be narrower and have poorer coverage when measurement error is small
#' est_cutt(ages=ages, sd=500, K=25000, alpha=0.05, method="mle") 
#' 
#' # Now repeat but with larger measurement error sd
#' ages5 = rcutt(20, 10000, K=25000, sd=5000)
#' 
#' # note this will run faster because it will call method-for a point estimate together with a 95% CI
#' est_cutt(ages=ages5, sd=5000, K=25000, alpha=0.05) 

est_cutt = function(ages, sd, K, df=NULL, alpha=0.05, q=c(lo=alpha/2,point=0.5,hi=1-alpha/2), method=if(mean(sd)/(K-min(ages))<0.1) "reginv" else "mle", paramInits=NULL, iterMax=1000)
{  

  if(method=="reginv")
  {
    # set up result list.
    theta = rep(NA,length(q))
    if(is.null(names(q)))
      names(theta)=paste0("q=",q)
    else
      names(theta)=names(q)

    result=list(theta=theta,q=q,error=theta,iter=theta,converged=theta)
    for(iQ in 1:length(q))
    {
      if(q[iQ]==0.5) regMethod="lm" else regMethod="rq"
      ftI = reginv_cutt(ages,sd,K,df,q=q[iQ],paramInits=paramInits,iterMax=iterMax,method=regMethod)
      result$theta[iQ] = ftI$theta[1]
      result$error[iQ] = ftI$error[1]
      result$df        = ftI$df[1]
      result$iter[iQ]  = ftI$iter[1]
      result$converged[iQ] = ftI$converged[1]
    }
    result$method="reginv"
    class(result) = "reginv"
  }
  else
  {
    result = mle_cutt(ages, sd, K, df=NULL, q=q)
  }
  class(result) = "est_cutt"
  result$call   = match.call()
  result$data   = data.frame(ages=ages,sd=sd)
  return(result)
}
