#' Confidence Interval for a Parameter using Regression Inversion
#'
#' Estimates a confidence interval for a parameter using regression inversion, given a model to simulate data from (as a function of this parameter)
#' and a function to compute a test statistic, given a sample of data.
#' @param data Dataset we want to use for inference 
#' @param getT A function to compute the statistic that will be used for inversion, as a function of \code{data}, the target parameter \code{theta}, and any additional arguments listed in \code{...}.
#' @param simulateData A function to simulate data from the assumed model, as a function of the target parameter \code{theta} and any additional arguments listed in \code{...}.
#' @param paramInits A vector of initial values of \code{theta} to try. Must have at least two unique values in it, but this algorithm is more stable if you give it a dozen or so.
#' @param q A scalar value between 0 and 1 specifying the quantile at which we want to solve for \code{theta}. Defaults to \code{0.5}
#' @param iterMax Maximum number of values of \code{theta} to generate new test statistics at.
#' @param eps Convergence tolerance for \code{theta}, defaults to \code{1.e-5}.
#' @param method The regression method employed for inversion. \code{'rq'} (default) uses quantile linear regression, \code{'rq2'} uses quantile quadratic regression,
#' \code{'wrq'} downweights influential values before using quantile linear regression, \code{'prob'} uses probit linear regression on an indicator variable for whether or not the observed test statistic \code{getT(data)} has been exceeded.
#' @param a A scalar value determining where to simulate the next value of \code{theta} at, relative to our best current estimate of it. \code{a=0} (default) simulates the next statistic at our current best estimate of \code{theta},
#' \code{a=1} simulates at the value the mean of all previous estimates equals our best estimate, intermediate values take a weighted mean of these two options. 
#' @param stats A data frame of previous results to be read in to assist in calculations. Should be used with care. Data frame should contain the test statistics (labelled as \code{T}) and the values of \code{theta} at which data were simulated. 
#' @param ... Further arguments passed through to \emph{both} the \code{getT} and \code{simulateData} functions. 
#' @details 
#' How well this works depends how well-chosen \code{getT} is (a maximum likelihood estimator, if available, is a great choice).
#' If there are nuisance parameters, their estimates could be passed as additional arguments (\code{...}) but a better approach is to write a function that can estimate these from \code{data} at any given input value of \code{theta}.
#' Then this function is called inside \code{simFn}, which in effect makes the simulation model a function of \code{theta} only.
#' 
#' @return This function returns an object of class "reginv" with the following components:
#'
#'  \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}. (If \code{q} was not specified as input, this defaults to the lower limit for a \code{100(1-alpha)}\% confidence interval, a point estimate at \code{q=0.5} ("best estimate" of extinction time), and an upper limit for a \code{100(1-alpha)}\% confidence interval.)}
#'  \item{B}{ a vector containing the number of Monte Carlo samples used to estimate each \code{theta}.}
#'  \item{q}{ the quantile used in estimation.}
#'  \item{call }{ the function call}
#' @examples
#' # Find the lower 2.5% confidence bound on the variance for a
#' # mean-zero Gaussian sample, using mean of squared observations as statistic
#' 
#' dat = rnorm(20, mean=0, sd=1) #simulating some random data
#' meanSq = function(x,n){sum(x^2)/n}
#' simNorm = function(x,n){rnorm(n,mean=0,sd=sqrt(x))}
#' 
#' reginv(dat, getT=meanSq, simulateData=simNorm, q=0.025, 
#'        paramInits=seq(0.5,2,length=20), n=length(dat))
#'        
#' # Compare to the exact value:
#' print( sum(dat^2)/qchisq(0.975,20) )

#' @export
#' @import stats
#' @importFrom quantreg rq
reginv = function(data, getT, simulateData, q=0.5, paramInits, iterMax=1000, eps=1.e-5, method="rq", a=0, stats = NULL, ...)
{
  t_obs = getT(data, ...)
  
  # get initial stats
  if(is.null(stats))
  {
    stats=data.frame(theta=NULL,T=NULL,thetaEst=NULL)
    for(i in 1:length(paramInits) )
    {
      newDat = simulateData(paramInits[i], ...)
      stats = rbind(stats, c(theta=paramInits[i], T = getT(newDat, ...), thetaEst=paramInits[i]) )
    }      
    names(stats)=c("theta","T","thetaEst")
  }
  
  # define getPred and getErr functions (to avoid repeated if statements in estimation)
  if(method=="rq"||method=="wrq"||method=="rq2")
  {
    q = 1-q # flipping quantile is needed when using this approach
    getPred = function(newtheta,qfit,t_obs,q)
    {
      return(predict(qfit,newdata=list(theta=newtheta))-t_obs)
    }
    getErr = function(newtheta,qfit,t_obs,q) # getting the SE of predictions 
    {
      prEst = try(predict(qfit,newdata=list(theta=newtheta),interval="confidence"),silent=TRUE)
      if(inherits(prEst,"try-error"))
        err = Inf
      else
      {
        if(coef(qfit)[2]>0) # find the value matching to observed t if decent fit
          err = (prEst[3]-prEst[2])/qnorm(0.975)/2/coef(qfit)[2]
        else
          err = Inf
      }
      if(is.na(err)) err=Inf
      return(err)
    }
  }
  else
  {
    getPred = function(newtheta,qfit,t_obs,q)
    {
      return(predict(qfit,newdata=list(theta=newtheta),type="response")-q)
    }
    getErr = function(newtheta,qfit,t_obs,q)
    {
      if(coef(qfit)[2]>0) # find the value matching to observed t if decent fit
      {
        prEst = predict(qfit,newdata=list(theta=newtheta),se.fit=TRUE)
        err = prEst$se.fit/coef(qfit)[2]
      }
      else 
        err = Inf
      return(err)
    }
  }
  
  iter = dim(stats)[1]
  res = updateTheta(stats, t_obs, q, method=method, getPred=getPred)
  isConverged = FALSE
  while(isConverged == FALSE & iter<iterMax)
  {
    iter=iter+1
    thetaNew = getThetaSim(iter,thetaEst=res$theta,thetaSims=stats$theta,a=a)
    newDat = simulateData(thetaNew, ...)
    Titer = getT(newDat, ...)
    if(is.na(Titer)|Titer==Inf) 
    {
      warning("getT not a number at current value of theta, resetting to median")
      thetaNew = median(stats$theta)
      newDat = simulateData(thetaNew, ...)
      Titer = getT(newDat, ...)
    }
    stats = rbind(stats, c(theta=thetaNew, T=Titer, thetaEst=res$theta) )
    res = updateTheta(stats, t_obs, q, qfit=res$qfit, method=method, getPred=getPred)
    err = getErr(res$theta,res$qfit,t_obs,q) / abs(res$theta)
    isConverged = err < eps
  }
  result = list(theta=res$theta,error=err,iter=iter,converged=isConverged,stats=stats,fit=res$qfit,call=match.call())
  names(result$theta)=paste0("q=",q)
  class(result)="reginv"
  return(result)
}

updateTheta = function(stats,t_obs,q,qfit=NULL,method="rq",getPred=getPred,screenStats=NULL)
  # dat is a dataframe containing theta and T (parameter and statistic)
{
  if(method=="wrq")
  {
    stats$wt = NA
    lm_ft = lm(T~theta,data=stats)
    infl = influence(lm_ft)$hat
    stats$wt[is.na(stats$T)==FALSE] = min(infl,2/length(stats$theta)) / infl
    wt=stats$wt # to avoid weird CRAN error
  }
  
  if(is.null(qfit))
  {
    ft = switch(method,
                "rq" = quantreg::rq(T~theta,tau=q,data=stats), #consider varying rq method, "fn" or "pfn"
                "rq2" = quantreg::rq(T~poly(theta,2,raw=TRUE),tau=q,data=stats), #consider varying rq method, "fn" or "pfn"
                #           "qgam" = qgam::qgam(T~s(theta), qu=q, data=stats),
                "wrq" = quantreg::rq(T~theta,tau=q,weights=wt,data=stats),
                glm(I(T>t_obs)~theta,family=binomial("probit"),data=stats)
    )
  }
  else
  {
    #update fit only using data points that are "local" to observed stat (predicted T within factor of screenStats of observed T)
    if(is.null(screenStats)==FALSE) 
    {
      prRatio = abs( log(predict(qfit,newdata=stats)/t_obs) )
      prRatio[is.na(prRatio)] = Inf
      OKstats = prRatio<log(screenStats)
      if(sum(OKstats)<10) #if this does not keep enough, take 5 closest
        OKstats=sort(prRatio,index.return=TRUE)$ix[1:10]
      statsWorking = stats[OKstats,]
    }
    else
      statsWorking=stats
    #    if(method=="qgam")   
    #      ft = update(qfit,data=stats,family=qfit$family)
    #    else
    ft = update(qfit,data=statsWorking)
  }
  if(method=="rq2")
  {
    # solve quadratic equation for larger solution if it exists and if stationary point is below first quartile
    Delta = coef(ft)[2]^2 - 4 * coef(ft)[3] * ( coef(ft)[1]-t_obs )
    #    statPt = -coef(ft)[2]/2/coef(ft)[3]
    #    cond = Delta>0 & statPt<quantile(stats$theta,0.25)
    cond = Delta>0
    if(cond)
      theta = ( -coef(ft)[2] + sqrt(Delta) ) / 2 / coef(ft)[3]
  }
  else
  {
    # solve linear equation if increasing slope
    if(method=="prob") t_obs=qnorm(q) # if using probit regression
    cond = coef(ft)[2]>0
    if(cond) # find the value matching to observed t if increasing fit
      theta = (t_obs-coef(ft)[1])/coef(ft)[2] 
  }
  if(cond == FALSE) # if fit is screwy use current best estimate + jitter until something changes!
    theta = stats$thetaEst[nrow(stats)]*(1+runif(1,-0.01,0.01))
  # don't accept any crazy shit... move them back towards the data
  # theta = noOutliers(theta,stats$theta)
  return(list(theta=theta,qfit=ft))
}

getThetaSim = function(iter,thetaEst,thetaSims,a)
{
  # a tells us how much to weight mean of all obs compared to current obs
  thetaNew = (a*iter+1-a)*thetaEst - a*sum(thetaSims)
  # put bounds on thetaNew so no crazy shit
  #  if(a>0) 
    thetaNext = noOutliers(thetaNew,thetaSims)
  # thetaNext=thetaNew
  return(thetaNext)
}
noOutliers = function(thetaNew,thetas)
# don't accept any crazy shit from theta - truncate using an asymmetric 3IQR rule
{
  quants = quantile(thetas,c(0.25,0.5,0.75))
  outRange = diff(quants)*3
  if ( thetaNew>quants[3]+outRange[2] )
    thetaNew = quants[3]+outRange[2]
  if ( thetaNew<quants[1]-outRange[1] )
    thetaNew = quants[1]-outRange[1]
  return(thetaNew)
}