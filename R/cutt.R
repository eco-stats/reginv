#' The Compound Uniform-Truncated T (CUTT) Distribution for Fossil Dates
#'
#' Density, distribution function, quantile function and random generation for the Compound Uniform-Truncated T (CUTT) distribution
#' for fossil dates with extinction time \code{theta}, oldest observable data \code{K} and measurement error \code{sd}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n Number of observations (numeric)
#' @param theta Time of extinction (numeric)
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param sd Measurement error standard deviations for fossils. Can be a vector. If not of length \code{n}, it is cycled through repeatedly until length \code{n}.
#' @param df numeric; degrees of freedom for the t-distribution used to model measurement error. Must be at least 2. Default (NULL) uses a Gaussian distribution.
#' @param log logical; if TRUE, densities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities \code{p} are \eqn{P[X\leq x]}{P(X<=x)} otherwise \eqn{P[X>x]}{P(X>x)}.
#' @param tol Numerical tolerance (defaults to \code{sqrt(.Machine$double.eps)}
#' @param nIter The maximum number of iterations to use estimating random numbers before resorting to \code{uniroot} (which is slower).
#' @param pMinus An amount to subtract from the distribution function in \code{pcutt} (used for root-finding in \code{rcutt} and \code{qcutt}).
#'
#' @details 
#' Given an extinction time \code{theta}, an upper limit \code{K} for the possible age of a fossil
#' that could be included in this dataset, and measurement error standard deviations \code{sd}, our procedure then assumes that:
#' \itemize{
#' \item For a fossil of a given age, measurement error estimating its age is normally distributed with mean zero and standard deviation as provided
#' \item That fossil dates are uniformly distributed over the interval of allowable dates
#' \item All observed fossil dates are less than \code{K}
#' }
#' This leads to a distribution function of the form
#' \deqn{C\left[(w-\theta)\Phi\left(\frac{w-\theta}{\sigma}\right)+\sigma^2\phi\left(\frac{w-\theta}{\sigma}\right)\right]}
#' where \eqn{C^{-1} = (K-\theta)\Phi\left(\frac{K-\theta}{\sigma}\right)+\sigma^2\phi\left(\frac{K-\theta}{\sigma}\right)}. 
#' 
#' @return \code{dcutt} gives the density, \code{pcutt} gives the distribution function, \code{qnorm} 
#' gives the quantile function, and \code{rcutt} generates random fossil dates.
#' 
#' The length of the result is determined by \code{n} for \code{rcutt}, and in other cases, by the length of the first argument.
#' \code{sd} can also be a vector but if its not length does not match that expected by the first argument it is cycled through until the desired
#' length is reached, with a warning.

##' @rdname cutt
##' @export
#' @examples
#' ages = rcutt(20, 10000, 25000, 1000) # simulate some random data
#' pcutt(ages, 10000, 25000, 1000) # find CDF of each value (approx uniform)
#' qcutt(c(0.25,0.75), 10000, 25000, 1000) # find first and third quartiles
#'
#' # plot the density function
#' w = seq(5000,26000,length=1000)
#' plot(w,dcutt(w,10000,25000,1000),type="l",ylab="pdf(w)")
#' # compare to density if measurement error came from t(4) distribution
#' plot(w,dcutt(w,10000,25000,1000,df=4),type="l",ylab="pdf(w,df=4)")
#' @aliases cutt rcutt dcutt pcutt qcutt
rcutt = function(n, theta, K, sd, df=NULL, tol=sqrt(.Machine$double.eps), nIter=50)
{
  if(length(n)>1) stop("'n' must be scalar")
  sdVec = checkParams(n, theta, K, df, sd)
  w = qcutt(p=runif(n), theta=theta, K=K, sd=sdVec, df=df, tol=tol, nIter=nIter)
  return(w)
}

##' @rdname cutt
##' @export
qcutt = function(p, theta, K, sd, df=NULL, tol=sqrt(.Machine$double.eps), nIter=50)
{
  # sort out the dimensions of inputs
  nP = length(p)
  sdVec = checkParams(nP,theta,K,df,sd)

  # get functions to compute F(e) and int e(f(e)) for measurement error e
  funs = getDFs(df,sd)
  
  # compute C and p/C
  Ksd  = (K-theta)/sdVec
  F.K  = funs$CDF(Ksd,df)
  f.K  = funs$fe(Ksd,df)
  Cinv = (K-theta) * F.K + sdVec * f.K
  pOnC = p * Cinv
  
  # now get cracking finding w
  q      = qOld = theta + p*(K-theta) # starting estimate (solution when sd=0)
  F.q    = f.q = qSD = rep(0,nP) 
  iter   = 0
  isDiff = sdVec!=0 #only do the below when sd is non-zero

  # check for silly p's, set to NaN and remove from consideration
  wrongP = p<0 | p>1
  if(any(wrongP))
  {
    q[wrongP]=NaN
    warning("NaNs produced")    
    isDiff[wrongP]=FALSE
  }
  while(any(isDiff) & iter<nIter)
  {
    qSD[isDiff]  = (qOld[isDiff]-theta)/sdVec[isDiff]
    F.q[isDiff]  = funs$CDF(qSD[isDiff],df)
    f.q[isDiff]  = funs$fe(qSD[isDiff],df)
    q[isDiff]    = theta + ( pOnC[isDiff] - sdVec[isDiff] * f.q[isDiff] ) / F.q[isDiff]
    q[isDiff]    = pmin(q[isDiff],K) #ensure it is less than K
    qDiff        = abs(q-qOld)
    isDiff       = qDiff>tol
    qOld[isDiff] = q[isDiff]
    iter         = iter+1
  }
  
  # in the occasional odd case this might not converge, in which case try using uniroot on cdf function pcutt
  if(iter==nIter)
  {
    for(iObs in which(isDiff))
    {
      eCrit = ifelse(is.null(df), qnorm(p[iObs]/2), qt(p[iObs]/2,df) )
      qTry = try( uniroot( pcutt, interval=c(theta+sdVec[iObs]*eCrit,K),tol=tol,theta=theta,sd=sdVec[iObs],K=K,df=df,pMinus=p[iObs],extendInt="upX") )

            if(inherits(qTry,"try-error")==FALSE)
      {
        q[iObs] = qTry$root
        isDiff[iObs] = FALSE
      }
    }
    if(any(isDiff)) warning(paste0("non-convergence for ",sum(isDiff)," observations"))
  }
  return(q)
}

##' @rdname cutt
##' @export
pcutt = function(q,theta,K,sd,df=NULL,lower.tail=TRUE,pMinus=0)
# function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # sort out the dimensions of inputs
  sdVec = checkParams(length(q),theta,K,df,sd)

  # get functions to compute F(e) and int e(f(e)) for measurement error e
  funs = getDFs(df,sd)
  
  # get CDF denominator
  Ksd  = (K-theta)/sdVec
  F.K  = funs$CDF(Ksd,df)
  f.K  = funs$fe(Ksd,df)
  Cinv = (K-theta) * F.K + sdVec * f.K
  # get CDF-u
  qSD  = (q-theta)/sdVec
  F.q  = funs$CDF(qSD,df)
  f.q  = funs$fe(qSD,df)
  cdf  = ( (q-theta) * F.q + sdVec * f.q ) / Cinv - pMinus
  cdf[q>K] = 1 - pMinus #set CDF to one if q>K
  if (lower.tail==FALSE) cdf=1-cdf
  return(cdf)
}

##' @rdname cutt
##' @export
dcutt = function(x,theta,K,sd,df=NULL,log=FALSE)
# function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # sort out the dimensions of inputs
  sdVec = checkParams(length(x),theta,K,df,sd)
  
  # get functions to compute F(e) and int e(f(e)) for measurement error e
  funs = getDFs(df,sd)

  # get CDF denominator
  Ksd  = (K-theta)/sdVec
  F.K  = funs$CDF(Ksd,df)
  f.K  = funs$fe(Ksd,df)
  Cinv = (K-theta) * F.K + sdVec * f.K
  
  # get CDF-u
  xSD  = (x-theta)/sdVec
  if(log)
  {
    pdf = funs$CDF(xSD,df,log.p=TRUE) - log(Cinv)
    pdf[x>K] = -Inf
  }
  else
  {
    pdf = funs$CDF(xSD,df)/Cinv
    pdf[x>K] = 0
  }
  return(pdf)
}

# checking if parameters are all the right dimensions
checkParams = function(n,theta,K,df,sd)
{
  if(length(theta)>1) stop("'theta' must be scalar")
  if(length(K)>1) stop("'K' must be scalar")
  if(is.null(df)==FALSE)
  {
    if(length(df)>1) stop("'df' must be scalar")
    if(df<2) stop("'df' must be at least 2")
  }
  #ensure sds is the right length
  sdVec=sd #return this if sd has length n
  nSD = length(sd)
  if(nSD==1) sdVec = rep(sd,n)
  if(nSD!=n & nSD!=1)
  {
    sdVec = rep(sd,length=n)
    warning("length of 'sd' is not equal to length of first argument so will extend 'sd' as needed.")
  }
  return(sdVec)
}

# define measurement error distribution function F(e) and integral of ef(e), used to compute cutt CDF
getDFs = function(df,sd)
{
  if(all(sd==0)) # if all sds are zero then returns 1 for CDF and 0 for pdf
  {
    CDFx = function(xSD,df,log.p=FALSE){ if(log.p) 0 else 1 }
    fex = function(xSD,df){ 0 }
  }
  else
  {
    if(is.null(df))
    {
      CDFx  = function(xSD,df,log.p=FALSE){pnorm(xSD, mean = 0, sd = 1, log.p=log.p)}
      fex  = function(xSD,df){dnorm(xSD, mean = 0, sd = 1)}
    }
    else
    {
      CDFx  = function(xSD,df,log.p=FALSE){ pt( xSD, df, log.p=log.p) }
      fex  = function(xSD,df){ dt(xSD,df) * df/(df-1) * (1+xSD^2/df) }
    }
  }
  return(list(CDF=CDFx,fe=fex))
}  
