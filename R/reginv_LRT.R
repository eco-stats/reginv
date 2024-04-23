#' Confidence Interval for Extinction Time using Bootstrap Inversion of a Likelihood Ratio Statistic
#'
#' Estimates a confidence interval for extinction time using bootstrap inversion, which finds a range of values for extinction time that
#' are plausible considering the value of the likelihood ratio statistic. This is done assuming fossil dates come from a
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
#' @param signroot logical, if set to \code{TRUE} (default) Will use a signed root likelihood ratio statistic, hence estimate asymmetric quantiles. If \code{FALSE} will use the likelihood ratio statistic, assuming symmetric quantiles.
#' @param iterMax Maximum number of simulated datasets to use to estimate extinction time (default 1000).
#' @param dfMin The minimum allowable value of the degrees of freedom, default=4.
#'
#' @details 
#' Given a vector of fossil ages \code{ages} and corresponding measurement error standard deviations \code{sd}, and an upper limit \code{K} for the possible age of a fossil
#' that could be included in this dataset, our procedure then assumes:
#' \itemize{
#' \item That, for a fossil of a given age, measurement error estimating its age follows a distribution that is Student's T multiplied by the provided standard deviation, truncated such that observed age is less than \code{K}
#' \item That, given a value for measurement error, fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until the value such that observed age is \code{K})
#' }
#' We then estimate extinction time by inversion of the (signed root) likelihood ratio statistic (LRT), that is, we find the estimate of extinction time \eqn{\theta}{t}
#' such that if this were the true extinction time then the chance of seeing a LRT less than the one obtained from sample \code{ages} is equal to \code{q}. The probability
#' is estimated via parametric bootstrap using parameters as estimated by \code{\link{mle_cutt}}. If \code{alpha} is specified this function will return three values:
#' the lower limit of the \code{100*(1-alpha)}\% confidence interval (solving at \code{q=alpha/2}), a point estimator for extinction time (solving at \code{q=0.5}), and an upper limit for the
#' confidence interval (solving at \code{q=1-alpha/2}). If a vector \code{q} is specified as input then the function solves for this vector instead. 
#' 
#' Note that because we are using simulation to estimate parameters and their confidence intervals, we will get slightly different answers on each run. 
#' 
#' It is assumed that \code{ages} has been specified with smaller values representing more recent specimens, for example, \code{ages} could be specified in years before present.
#' If there is interest in estimating speciation or invasion time, data would only need to be reordered so that smaller values represent older specimens (e.g. by multiplying by negative one).
#' 
#' @return This function returns an object of class \code{est_cutt} with the following components:
#'
#'  \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}. (If \code{q} was not specified as input, this defaults to a bias-corrected maximum likelihood estimate \code{point}, and together with confidence interval limits \code{lo} and \code{hi}.)}
#'  \item{se}{ the estimated standard error of the MLE.}
#'  \item{df}{ the estimated degrees of freedom of the Student's T distribution for measurement error.}
#'  \item{data}{ a data frame with the data used in analysis, in columns labelled as \code{ages} and \code{sd}.}
#'  \item{method}{ \code{bootlrt}.}
#'  \item{call }{ the function call.}
#' @seealso est_cutt, mle_cutt, cutt
#' @export
#' @examples
#' ages = rcutt(20, 10000, K=25000, sd=500) #simulating some random data
#' 
#' # for a point estimate together with a 95% CI (only 200 iterations used so it runs quickly)
#' bootlrt_cutt(ages=ages, sd=500, K=25000, alpha=0.05, iterMax=200) 
#' 
#' # compare to estimates using asymptotic likelihood inference, which tend to
#' # be narrower and have poorer coverage (they miss the true value too often
#' # when n or sd is small):
#' mle_cutt(ages=ages, sd=500, K=25000, alpha=0.05) 
#' @export
bootlrt_cutt = function(ages, sd, K, df=NULL, alpha=0.05, q=c(lo=alpha/2,point=0.5,hi=1-alpha/2), signroot=TRUE, iterMax=1000, dfMin=4)
{  
  # set up result list.
  theta = rep(NA,length(q))
  if(is.null(names(q)))
    names(theta)=paste0("q=",q)
  else
    names(theta)=names(q)

  # get MLE and define get_LRTi function
  if(is.null(df))
  {
    mles = getJointMLE(ages=ages, theta=min(ages), sd=sd, K=K, dfMin=dfMin)
    thetaMLE = list(par=mles$par[1],value=mles$value)
    dfOut = 1/mles$par[2]
    vr = try( solve(-mles$hessian) )
    if(inherits(vr,"try-error")) vr=matrix(NaN,2,2)
    SE = if(is.nan(vr[1,1])) 0 else as.vector(sqrt(vr[1,1]))
    
    # define get_LRTi function
    get_LRTi = function(theta0, ages, K, sd, df=dfOut, const=0, dfMin=dfMin)
    {
      ll1  = getJointMLE(ages=ages, theta=theta0, sd=sd, K=K, df=df)
      ll0  = getDF( ages, theta=theta0, sd=sd, K=K, dfInvInit=1/df, dfMin=dfMin )
      return(sign(theta0-ll1$par[1])*sqrt(-2*(ll0$value-ll1$value)) - const)
    }
  }
  else
  {
    dfOut=df
    thetaMLE = getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=df)
    SE = as.vector(1/sqrt(-thetaMLE$hessian))
    
    # define get_LRTi function
    get_LRTi = function(theta0, ages, K, sd, df=dfOut, const=0, dfMin=dfMin)
      {
      ll1  = getThetaMLE(ages, theta=theta0, sd=sd, K=K, df=df )
      ll0  = cutt_LogLik(theta0,ages,sd=sd,K=K,df=df)
      return(sign(theta0-ll1$par)*sqrt(-2*(ll0$value-ll1$value)) - const)
    }
  }

  # get LRs
  n = length(ages)
  thMLE = thetaMLE$par
  # resample iterMax times
#  LRs = rep(NA,iterMax)
#  for(b in 1:iterMax)
#  {
#    ageStar = rcutt(n,theta=thMLE,K=K,sd=sd,df=dfOut)
#    LRs[b]  = get_LRTi(thMLE,ages=ageStar,K=K,sd=sd,df=dfOut, dfMin=dfMin)
#  }
  # get sample quantiles of LRT
#  qLR = quantile(LRs, q, na.rm=TRUE)
  thetaWorking = mle_cutt(ages=ages, sd=sd, K=K, df=df, q=q)
  thetaWorking$theta = pmin(thetaWorking$theta,K) #don't let max go crazy
    
  dir="upX" #increasing function of theta

  # set search limits so that we look above MLE if q>0.5 and below otherwise 
  is_SE_bad   = is.nan(SE) | is.infinite(SE) | SE==0
  searchLim   = ifelse( is_SE_bad, IQR(ages)*0.5, SE*5 )
  qLo = qHi   = rep(thMLE,length(q))
  qLo[q<=0.5] = thMLE-searchLim
  qHi[q>=0.5] = min(thMLE,K)+searchLim
  
  q2 = 1-2*pmin(q,1-q) #two-sided quantile

  # get results for each value of q
  for (iQ in 1:length(q))
  {
#   to check problematic cases
#    ths = seq(qLo[2],qHi[2],length=100)
#    lrs=rep(NA,100)
#    for(i in 1:100) lrs[i] = get_LRTi(ths[i],ages=ages,K=K,sd=sd,dfOut,const=0,dfMin)
#    plot(lrs~ths)

    if(q[iQ]==0.5)
      theta[iQ] = thMLE
    else
    {
      
      # get estimates to use in simulation
      dfWorking  = getDF( ages, theta=thetaWorking$theta[iQ], sd=sd, K=K, dfInvInit=1/dfOut, dfMin=dfMin )$par
      
      LRs = rep(NA,iterMax)
      for(b in 1:iterMax)
      {
        ageStar = rcutt(n,theta=thetaWorking$theta[iQ],K=K,sd=sd,df=dfWorking)
        LRs[b]  = get_LRTi(thetaWorking$theta[iQ],ages=ageStar,K=K,sd=sd,df=dfWorking, dfMin=dfMin)
      }
      # get sample quantiles of LRT

      qLR = ifelse(signroot,
                   quantile(LRs, q[iQ], na.rm=TRUE),      
                   quantile(abs(LRs), q2[iQ], na.rm=TRUE)*sign(q[iQ]-0.5)
      )
      thLim = try( uniroot(f=get_LRTi, interval=c(qLo[iQ],qHi[iQ]), ages=ages, K=K, sd=sd, df=dfWorking, const=qLR, dfMin=dfMin, extendInt=dir) )
      
#      thLim = try( uniroot(f=get_LRTi, interval=c(qLo[iQ],qHi[iQ]), ages=ages, K=K, sd=sd, df=dfWorking, const=qLR[iQ], dfMin=dfMin, extendInt=dir) )
      if(inherits(thLim,"try-error"))
        theta[iQ] = thetaWorking$theta[iQ]
      else
        theta[iQ] = thLim$root
    }
  }
  result = list(theta = theta, se=thetaMLE$se, df=dfOut, data = data.frame(ages=ages, sd=sd), method = "bootlrt", call = match.call() )
  class(result) = "est_cutt"
  return(result)
}


#' @export
reginv_LRT = function(ages, sd, K, df=NULL, alpha=0.05, q=c(lo=alpha/2,point=0.5,hi=1-alpha/2), iterMax=1000, method="rq")
{  
  
  # set up result list.
  theta = rep(NA,length(q))
  if(is.null(names(q)))
    names(theta)=paste0("q=",q)
  else
    names(theta)=names(q)
  result=list(theta=theta,q=q,error=theta,iter=theta,converged=theta)
  

  if(is.null(df)==FALSE) dfi=df

  # get results for each value of q
  for (iQ in 1:length(q))
  {
    ft.mle = mle_cutt(ages=ages, sd=sd, K=K, df=df, q=q[iQ])
    if(q[iQ]==0.5)
    {
      result$theta[[iQ]] = ft.mle$theta[1]
      result$se = ft.mle$se
      result$df = ft.mle$df
    }
    else
    {
      # get paramInits
      is_SE_bad = is.nan(ft.mle$se) | is.infinite(ft.mle$se) | ft.mle$se==0
      stepSize = ifelse( is_SE_bad, IQR(ages)*0.1, ft.mle$se )
      if(ft.mle$theta[1]+5*stepSize>K) #make sure paramInits don't exceed K
      {
        if(ft.mle$theta[1]<K)
          paramInits = seq( ft.mle$theta[1] -5* stepSize, K, length=100) 
        else
          paramInits = seq( mle_cutt(ages=ages, sd=sd, K=K, df=df, alpha=NULL)$theta-5*stepSize,K,length=100 )
      }
      else
        paramInits = ft.mle$theta[1] + stepSize*seq(-5,5,length=100) 
      
      n <- length(ages)
      stats=data.frame(theta=NULL,Tnew=NULL,Tobs=NULL, df=NULL)
      for(i in 1:length(paramInits) )
      {
        if( is.null(df) )  dfi = getDF(ages, theta=paramInits[i], sd=sd, K=K)$par
        statObs=get_LRT(paramInits[i],ages, K=K, sd=sd, df=dfi)
        newDat = simFn_cutt(theta=paramInits[i],K=K,sd=sd,df=dfi,n=n,dat=ages)
        if( is.null(df) )  dfn = getDF(newDat, theta=paramInits[i], sd=sd, K=K)$par
        statNew=get_LRT(paramInits[i],newDat, K=K, sd=sd, df=dfn)
        stats = rbind(stats, c(theta=paramInits[i], Tnew = statNew, Tobs=statObs, df=dfn) )
      }      
      names(stats)=c("theta","Tnew","Tobs","df")
      
      # now iterate to convergence
      iter = dim(stats)[1]
      res = updateTheta_LRT(stats, q[iQ], K=K, sd=sd)
      while(iter<iterMax)
      {
        iter=iter+1
        
        if( is.null(df) )  dfi = getDF(ages, theta=res$theta, sd=sd, K=K)$par
        statObs=get_LRT(res$theta,ages, K=K, sd=sd, df=dfi)
        newDat = simFn_cutt(theta=res$theta,K=K,sd=sd,df=dfi,n=n,dat=ages)
        if( is.null(df) )  dfn = getDF(newDat, theta=res$theta, sd=sd, K=K)$par
        statNew=get_LRT(res$theta,newDat, K=K, sd=sd, df=dfn)
        stats = rbind(stats, c(theta=res$theta, Tnew = statNew, Tobs=statObs, df=dfn) )
        
        res = updateTheta_LRT(stats, q[iQ], K=K, sd=sd, qfit=res$qfit)
      }
      result$theta[[iQ]] = res$theta
    }
  }
  result$data = data.frame(ages=ages, sd=sd)
  result$method="reginv"
  result$call <- match.call()
  result$stats=stats
  class(result)="reginv"
  class(result)="est_cutt"
  return(result)
}

get_LRT = function(theta0, ages, K, sd, df)
{
  thetaMLE = getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=df)
  ll0 = cutt_LogLik(theta0,ages,sd=sd,K=K,df=df)
  LRstat = sign(theta0-thetaMLE$par) * sqrt( 2*(thetaMLE$value-ll0) )
}

updateTheta_LRT = function(stats,q,K,sd,qfit=NULL)
  # dat is a dataframe containing theta and T (parameter and statistic)
{
  if(is.null(qfit))
    ft = quantreg::rq(Tnew~theta,tau=q,data=stats) #consider varying rq method, "fn" or "pfn"
  else
  {
    statsWorking = stats
    ft = update(qfit,data=statsWorking)
  }
  # find difference in (predicted) statistics and look for a zero
  Tdiff   = stats$Tobs - coef(ft)[1] - stats$theta*coef(ft)[2]
  posZero = min(Tdiff[Tdiff>=0])
  negZero = max(Tdiff[Tdiff<=0])
  
  whichPos = which(Tdiff==posZero)
  whichNeg = which(Tdiff==negZero)
  if(length(whichPos)>1) whichPos = whichPos[which.min(stats$theta[whichPos])]
  if(length(whichNeg)>1) whichNeg = whichNeg[which.max(stats$theta[whichNeg])]
  thetaNeg = stats$theta[whichNeg]
  thetaPos = stats$theta[whichPos]
  TNeg     = Tdiff[whichNeg]
  TPos     = Tdiff[whichPos]
  
  # linearly interpolate nearest values to zero if there is a solution
  if(length(whichPos)==1 & length(whichNeg)==1)
  {
    theta    = thetaNeg - TNeg/(TPos-TNeg)*(thetaPos-thetaNeg)
    if(TPos-TNeg==0) theta=mean(c(thetaNeg,thetaPos)) # if both are zeros average them
  }

  #if no sign change then extend the range a bit...
  if(length(whichNeg)<1) theta=thetaPos-0.25*diff(range(stats$theta))
  if(length(whichPos)<1) theta=thetaNeg+0.25*diff(range(stats$theta))
  theta = min(theta,K+mean(sd))

  return(list(theta=theta,qfit=ft))
}
