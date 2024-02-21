#'@export
print.reginv <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n")
  print(x$theta)
  cat(" CI estimated using", x$method, "method\n")
}

#'@export
print.est_cutt <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  if(is.null(x$q)==FALSE)
  {
    cat("\n   CI:\n")
    print(x$theta)
  }
  else
    cat("\n  MLE:", x$theta[1],"\n")
  cat(" CI estimated using", x$method, "method\n")
}

#'@export
simulate.est_cutt = function (object, nsim = 1, seed = NULL, ...) 
{
  y=getCUTTparams(object)
  out = rcutt(nsim*length(y$ages),theta=y$theta,K=y$K,sd=rep(y$sd,nsim),df=y$df)
  if(nsim>1)
  {
    out = matrix(out,ncol=nsim)
    rownames(out)=names(y$ages)
    colnames(out)=paste("sim_",1:nsim)
  }
  else
    names(out)=names(y$ages)
  return(out)
}  


#' Quantile-Quantile Plots with Global Simulation Envelopes
#'
#' Produces a QQ plot from fitted object \code{y} with global envelopes constructed
#' by simulation. Global envelopes are constructed using the \code{GET} package for 
#' simultaneous control of error rates over the whole plot. For details see \code{\link[ecostats]{qqenvelope}}
#'
#' @param y is a fitted object.
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 199 datasets.
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param ylab \code{y} axis label (if a plot is produced).
#' @param ... further arguments sent through to \code{plot}.
#' @export
qqenvelope = function(y, n.sim=199, conf.level=0.95, ylab="Sample Quantiles", ...) UseMethod("qqenvelope")

#'@export
qqenvelope.default = function(y, n.sim=199, conf.level=0.95, ylab="Sample Quantiles", ...) ecostats::qqenvelope(y, n.sim=199, conf.level=0.95, ylab="Sample Quantiles", ...)

  
getCUTTparams = function(y)
# getting CUTT parameters out of a est_cutt object
{
  if(is.null(y$q))
  {
    theta = y$theta
  }
  else
    theta = y$theta[names(y$theta)=="point"]
  K = y$call[names(y$call)=="K"][[1]]
  return(list(ages=y$data$ages,theta=theta,K=K,sd=y$data$sd,df=y$df))
}