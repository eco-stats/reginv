% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootlrt_cutt.R
\name{bootlrt_cutt}
\alias{bootlrt_cutt}
\title{Confidence Interval for Extinction Time using Bootstrap Inversion of a Likelihood Ratio Statistic}
\usage{
bootlrt_cutt(
  ages,
  sd,
  K,
  df = NULL,
  alpha = 0.05,
  q = c(lo = alpha/2, point = 0.5, hi = 1 - alpha/2),
  signroot = TRUE,
  iterMax = 1000,
  dfMin = 4,
  ncpus = NULL
)
}
\arguments{
\item{ages}{Numeric vector of fossil ages, with smaller values being more recent fossil ages.}

\item{sd}{Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).}

\item{K}{Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{K} is
close to the age of the oldest fossil.}

\item{df}{Numeric; degrees of freedom for the t-distribution used to model measurement error. Must be greater than 2. Uses a Gaussian distribution if \code{df=Inf}.
Default (\code{NULL}) estimates \code{df} from the data.}

\item{alpha}{Numeric between 0 and 1. Used to find a 100(1-\code{alpha})\% confidence interval. Defaults to 0.05 (95\% confidence intervals)}

\item{q}{Numeric vector of values between 0 and 1, specifying the quantiles at which we want to solve for extinction time. Defaults to \code{c(alpha/2,0.5,1-alpha/2)},
which gives the limits of a 100(1-\code{alpha})\% confidence interval and a point estimate obtained by solving at 0.5. If \code{q} is specified it overrides any input for \code{alpha}.}

\item{signroot}{logical, if set to \code{TRUE} (default) Will use a signed root likelihood ratio statistic, hence estimate asymmetric quantiles. If \code{FALSE} will use the likelihood ratio statistic, assuming symmetric quantiles.}

\item{iterMax}{Maximum number of simulated datasets to use to estimate extinction time (default 1000).}

\item{dfMin}{The minimum allowable value of the degrees of freedom, default=4.}

\item{ncpus}{The number of cores used to do bootstrapping. Default \code{NULL} will use two less than the maximum available.}
}
\value{
This function returns an object of class \code{est_cutt} with the following components:

 \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}. (If \code{q} was not specified as input, this defaults to a bias-corrected maximum likelihood estimate \code{point}, and together with confidence interval limits \code{lo} and \code{hi}.)}
 \item{se}{ the estimated standard error of the MLE.}
 \item{df}{ the estimated degrees of freedom of the Student's T distribution for measurement error.}
 \item{data}{ a data frame with the data used in analysis, in columns labelled as \code{ages} and \code{sd}.}
 \item{method}{ \code{bootlrt}.}
 \item{call }{ the function call.}
}
\description{
Estimates a confidence interval for extinction time using bootstrap inversion, which finds a range of values for extinction time that
are plausible considering the value of the likelihood ratio statistic. This is done assuming fossil dates come from a
a compound uniform-truncated T distribution, which accounts for sampling error (the fact that the most recent fossil date is not necessarily
the most recent time that the species was extant) and measurement error (error dating fossils).
Usually takes a while to run (but less than a minute), because it is doing a lot of computation behind the scenes.
}
\details{
Given a vector of fossil ages \code{ages} and corresponding measurement error standard deviations \code{sd}, and an upper limit \code{K} for the possible age of a fossil
that could be included in this dataset, our procedure then assumes:
\itemize{
\item That, for a fossil of a given age, measurement error estimating its age follows a distribution that is Student's T multiplied by the provided standard deviation, truncated such that observed age is less than \code{K}
\item That, given a value for measurement error, fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until the value such that observed age is \code{K})
}
We then estimate extinction time by inversion of the (signed root) likelihood ratio statistic (LRT), that is, we find the estimate of extinction time \eqn{\theta}{t}
such that if this were the true extinction time then the chance of seeing a LRT less than the one obtained from sample \code{ages} is equal to \code{q}. The probability
is estimated via parametric bootstrap using parameters as estimated by \code{\link{mle_cutt}}. If \code{alpha} is specified this function will return three values:
the lower limit of the \code{100*(1-alpha)}\% confidence interval (solving at \code{q=alpha/2}), a point estimator for extinction time (solving at \code{q=0.5}), and an upper limit for the
confidence interval (solving at \code{q=1-alpha/2}). If a vector \code{q} is specified as input then the function solves for this vector instead. 

Note that because we are using simulation to estimate parameters and their confidence intervals, we will get slightly different answers on each run. This code is fairly computationally intensive but to shorten compute times
it makes use of multiple processors, as set using \code{ncpus}.

It is assumed that \code{ages} has been specified with smaller values representing more recent specimens, for example, \code{ages} could be specified in years before present.
If there is interest in estimating speciation or invasion time, data would only need to be reordered so that smaller values represent older specimens (e.g. by multiplying by negative one).
}
\examples{
ages = rcutt(20, 10000, K=25000, sd=500) #simulating some random data

# for a point estimate together with a 95\% CI (only 200 iterations used so it runs quickly)
bootlrt_cutt(ages=ages, sd=500, K=25000, alpha=0.05, iterMax=200, ncpus=2) 

# compare to estimates using asymptotic likelihood inference, which tend to
# be narrower and have poorer coverage (they miss the true value too often
# when n or sd is small):
mle_cutt(ages=ages, sd=500, K=25000, alpha=0.05) 
}
\seealso{
est_cutt, mle_cutt, cutt
}
