#' The Uniform-Gaussian Distribution for Fossil Dates
#'
#' Density, distribution function, quantile function and random generation for the Uniform-Gaussian distribution
#' for fossil dates with extinction time \code{theta}, oldest observable data \code{K} and measurement error \code{sd}.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n Number of observations (numeric)
#' @param theta Time of extinction (numeric)
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
#' close to the age of the oldest fossil.
#' @param sd Measurement error standard deviations for fossils. Can be a vector. If not of length \code{n}, it is cycled through repeatedly until length \code{n}.
#' @param log.p logical; if TRUE, densities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities \code{p} are \eqn{P[X\leq x]}{P(X<=x)} otherwise \eqn{P[X>x]}{P(X>x)}.
#' @param tol Numerical tolerance (defaults to \code{sqrt(.Machine$double.eps)}
#' @param nIter The maximum number of iterations to use estimating random numbers before resorting to \code{uniroot} (which is slower).
#' @param pMinus An amount to subtract from the distribution function in \code{pfossil} (used for root-finding in \code{rfossil} and \code{qfossil}).
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
#' @return \code{dfossil} gives the density, \code{pfossil} gives the distribution function, \code{qnorm} 
#' gives the quantile function, and \code{rfossil} generates random fossil dates.
#' 
#' The length of the result is determined by \code{n} for \code{rfossil}, and in other cases, by the length of the first argument.
#' \code{sd} can also be a vector but if its not length does not match that expected by the first argument it is cycled through until the desired
#' length is reached, with a warning.

##' @rdname fossil
##' @export
#' @examples
#' ages = rfossil(20, 10000, 25000, 1000) # simulate some random data
#' pfossil(ages, 10000, 25000, 1000) # find CDF of each value (approx uniform)
#' qfossil(c(0.25,0.75), 10000, 25000, 1000) # find first and third quartiles
#'
#' # plot the density function
#' w = seq(5000,26000,length=1000)
#' plot(w,dfossil(w,10000,25000,1000),type="l",ylab="pdf(w)")
#' @aliases fossil rfossil dfossil pfossil qfossil
rfossil = function(n, theta, K, sd, tol=sqrt(.Machine$double.eps), nIter=50)
{
  if(length(theta)>1) stop("'theta' must be scalar")
  if(length(K)>1) stop("'K' must be scalar")
  if(length(n)>1) stop("'n' must be scalar")
  #ensure sds is the right length
  nSD = length(sd)
  if(nSD!=n & nSD!=1)
  {
    sd = rep(sd,length=n)
    warning("length of 'sd' is not equal to 'n' so will extend 'sd' as needed.")
  }
  
  w = qfossil(p=runif(n), theta=theta, K=K, sd=sd, tol=tol, nIter=nIter)
  return(w)
}

##' @rdname fossil
##' @export
qfossil = function(p, theta, K, sd, tol=sqrt(.Machine$double.eps), nIter=50)
{
  # sort out the dimensions of inputs
  if(length(theta)>1) stop("theta must be scalar")
  if(length(K)>1) stop("K must be scalar")
  nSD = length(sd)
  nP  = length(p)
  if(nSD==1)
  {
    sd = rep(sd,length=nP)
    nSD = nP # to avoid below warning for scalar sd
  }
  if(nSD!=nP)
  {
    sd = rep(sd,length=nP)
    warning("length of 'sd' is not equal to length of 'p' so will extend 'sd' as needed.")
  }

  # compute C and p/C
  F.K  = pnorm(K - theta, mean = 0, sd = sd)
  f.K  = dnorm(K - theta, mean = 0, sd = sd)
  Cinv = (K-theta) * F.K + sd^2 * f.K
  pOnC = p * Cinv
  
  # now get cracking finding w
  q      = qOld = theta + p*(K-theta) # starting estimate (solution when sd=0)
  F.q    = f.q = rep(0,nP) 
  iter   = 0
  isDiff = sd!=0 #only do the below when sd is non-zero
  while(any(isDiff) & iter<nIter)
  {
    F.q[isDiff]  = pnorm(qOld[isDiff]-theta, mean = 0, sd = sd[isDiff])
    f.q[isDiff]  = dnorm(qOld[isDiff]-theta, mean = 0, sd = sd[isDiff])
    q[isDiff]    = theta + ( pOnC[isDiff] - sd[isDiff]^2 * f.q[isDiff] ) / F.q[isDiff]
    q[isDiff]    = pmin(q[isDiff],K) #ensure it is less than K
    qDiff        = abs(q-qOld)
    isDiff       = qDiff>tol
    qOld[isDiff] = q[isDiff]
    iter         = iter+1
  }
  
  # in the occasional odd case this might not converge, in which case try using uniroot on cdf function pfossil
  if(iter==nIter)
  {
    for(iObs in which(isDiff))
    {
      qTry = try( uniroot( pfossil, interval=c(qnorm(sqrt(tol),mean=0,sd=sd[iObs]),K-theta),tol=tol,n=1,theta=theta,sd=sd[iObs],K=K,pMinus=p[iObs],extendInt="upX") )
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

##' @rdname fossil
##' @export
pfossil = function(q,theta,K,sd,n=length(sd),lower.tail=TRUE,pMinus=0)
# function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # sort out the dimensions of inputs
  if(length(theta)>1) stop("theta must be scalar")
  if(length(K)>1) stop("K must be scalar")
  nSD = length(sd)
  nQ  = length(q)
  if(nSD!=nQ & nSD!=1)
  {
    sd = rep(sd,length=nQ)
    warning("length of 'sd' is not equal to length of 'q' so will extend 'sd' as needed.")
  }
  if (lower.tail==TRUE) q=1-q

  # get CDF denominator
  F.K  = pnorm(K - theta, mean = 0, sd = sd)
  f.K  = dnorm(K - theta, mean = 0, sd = sd)
  Cinv = (K-theta) * F.K + sd^2 * f.K
  
  # get CDF-u
  F.q  = pnorm(q-theta, mean = 0, sd = sd)
  f.q  = dnorm(q-theta, mean = 0, sd = sd)
  cdf  = ( (q-theta) * F.q + sd^2 * f.q ) / Cinv - pMinus
  cdf[q>K] = 1 - pMinus
  return(cdf)
}

##' @rdname fossil
##' @export
dfossil = function(x,theta,K,sd,log.p=FALSE)
  # function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # sort out the dimensions of inputs
  if(length(theta)>1) stop("theta must be scalar")
  if(length(K)>1) stop("K must be scalar")
  nSD = length(sd)
  nX  = length(x)
  if(nSD!=nX & nSD!=1)
  {
    sd = rep(sd,length=nX)
    warning("length of 'sd' is not equal to length of 'x' so will extend 'sd' as needed.")
  }

  # get CDF denominator
  F.K  = pnorm(K - theta, mean = 0, sd = sd)
  f.K  = dnorm(K - theta, mean = 0, sd = sd)
  Cinv = (K-theta) * F.K + sd^2 * f.K
  
  # get CDF-u
  if(log.p)
  {
    pdf = pnorm(x-theta, mean = 0, sd = sd, log.p=TRUE) - log(Cinv)
    pdf[x>K] = -Inf
  }
  else
  {
    pdf = pnorm(x-theta, mean = 0, sd = sd)/Cinv
    pdf[x>K] = 0
  }
  return(pdf)
}
