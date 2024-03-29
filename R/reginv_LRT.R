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
        if( is.null(df) )  dfi = reginv:::getDF(ages, theta=paramInits[i], sd=sd, K=K)$par
        statObs=get_LRT(paramInits[i],ages, K=K, sd=sd, df=dfi)
        newDat = reginv:::simFn_cutt(theta=paramInits[i],K=K,sd=sd,df=dfi,n=n,dat=ages)
        if( is.null(df) )  dfn = reginv:::getDF(newDat, theta=paramInits[i], sd=sd, K=K)$par
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
        
        if( is.null(df) )  dfi = reginv:::getDF(ages, theta=res$theta, sd=sd, K=K)$par
        statObs=get_LRT(res$theta,ages, K=K, sd=sd, df=dfi)
        newDat = reginv:::simFn_cutt(theta=res$theta,K=K,sd=sd,df=dfi,n=n,dat=ages)
        if( is.null(df) )  dfn = reginv:::getDF(newDat, theta=res$theta, sd=sd, K=K)$par
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

get_LRT = function(theta0,ages, K, sd, df)
{
  thetaMLE = reginv:::getThetaMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=df)
  ll0 = reginv:::cutt_LogLik(theta0,ages,sd=sd,K=K,df=df)
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
