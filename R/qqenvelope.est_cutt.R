#' Quantile-Quantile Plots with Global Simulation Envelopes for the CUTT distribution
#'
#' Produces a QQ plot from a \code{\link{est_cutt}} object that was fitted to a set of 
#' fossil ages, compared to "theoretical" quantiles, with global envelopes constructed
#' by simulation. Global envelopes are constructed using the \code{GET} package for 
#' simultaneous control of error rates over the whole plot.
#'
#' @param y is a \code{\link{est_cutt}} object.
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 199 datasets.
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param type the type of global envelope to construct, see 
#'  \code{\link[GET]{global_envelope_test}} for details. Default \code{"st"} uses 
#'  studentized envelope tests to protect for unequal variance, which has performed well
#'  in simulations in this context.
#' @param ylab \code{y} axis label (if a plot is produced).
#' @param main the plot title (if a plot is produced).
#' @param xlab \code{x} axis label (if a plot is produced).
#' @param ylab \code{y} axis label (if a plot is produced).
#' @param col color of points
#' @param line.col color of the lines on the diagnostic plot, defaults to "olivedrab", because it's cool.
#' @param envelope.col color of the global envelope around the expected trend. All data points should always stay within this envelope
#' (and will for a proportion \code{conf.level} of datasets satisfying model assumptions).
#' @param plot.it logical. Should the result be plotted? If not, a list of analysis outputs is returned, see \emph{Value}.
#' @param ... further arguments sent through to \code{plotenvelope} and \code{plot}.
#' 
#' @details
#' This function constructs a qqplot to check if a set of fossil ages are consistent with a \code{\link{cutt}}
#' distribution, whose parameters were estimated via \code{\link{est_cutt}}. A global simulation envelope
#' is included, estimated by simulating multiple realizations from the CUTT distribution.
#' All data points should lie if assumptions are satisfied, and will do so for a proportion \code{conf.level} of
#' datasets which satisfy their assumptions.
#' 
#' The envelope is global, constructed using the \code{\link[GET]{GET-package}}. So if \emph{any} data points lie outside the
#' envelope we have evidence that assumptions are not satisfied. 
#' The \code{\link[GET]{GET-package}} was originally constructed to place envelopes around functions, motivated by
#' the problem of diagnostic testing of spatial processes (Myllymäki et al 2017), but it can equally
#' well be applied here, by treating sorted residuals as point-wise evaluations of a function.
#' 
#' @return a qqplot with simulation envelope is returned, and additionally:
#' \item{x}{a vector of theoretical quantiles from the fitted CUTT distribution, sorted from
#'  smallest to largest}
#' \item{y}{a vector of observed fossil ages, sorted from smallest to largest}
#' \item{lo}{lower bounds on the global simulation envelope for fossil ages}
#' \item{hi}{upper bounds on the global simulation envelope for fossil ages}
#' \item{p.value}{A \emph{P}-value for the test that model assumptions are correct,
#'  using a 'parametric bootstrap' test, based on how far fossil ages depart from
#'  the values expected of them if CUTT model assumptions were satisfied.}
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' @references 
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017), Global envelope tests for spatial processes. J. R. Stat. Soc. B, 79: 381-404. doi:10.1111/rssb.12172
#' Warton DI (2022) Eco-Stats - Data Analysis in Ecology, from \emph{t}-tests to multivariate abundances. Springer, ISBN 978-3-030-88442-0
#' 
#' @seealso \code{\link{cutt}}, \code{\link{est_cutt}}, \code{\link{qqnorm}}, \code{\link{qqenvelope}} 
#' @examples
#' deer = c(14215, 14235, 14160, 14096, 14060, 13985, 13885, 13864, 13854, 13845, 13676, 13664,
#'  13624, 13644, 13590, 13576, 13491, 13445, 13431, 13400, 13379, 13383, 13358, 13349, 13347,
#'  13340, 13259, 13193, 13193, 13190, 13182, 13034, 13022, 12957, 12826, 12840, 12860, 12791, 
#'  12780, 12693, 12629, 12740, 12592, 12433)
#' deerSds = c(105.75, 78.5, 69.5, 102, 131.25, 120, 97.25, 89.5, 89.5, 95, 78.5, 74.75, 132.75, 
#'  94, 116.75, 119, 54.75, 119.5, 112.75, 60, 90.75, 58.25, 93, 92, 65.75, 102.25, 284.25, 
#'  102.25, 62.75, 75, 95.75, 70, 55.75, 114, 74.75, 110.25, 205.75, 90.75, 97.25, 66.75, 78.25, 
#'  307.75, 68, 96.25)
#' 
#' # get a maximum likelihood estimator and approx CIs
#' deerFt = est_cutt(ages=deer, sd=deerSds, K=14250, method="mle") 
#' 
#' # do a qq plot with simulation envelope around CUTT distribution fit to deer data:
#' qqenvelope(deerFt) 
#' @aliases qqenvelope.est_cutt
#' @import ecostats
#' @importFrom grDevices adjustcolor
#' @export
qqenvelope.est_cutt = function (y, n.sim=199, conf.level=0.95, ylab = "Sample quantiles", type="st", main = "Quantile Plot", 
                                xlab = "Theoretical Quantiles", col=NULL, line.col="olivedrab", 
                                envelope.col = adjustcolor(line.col, 0.1), plot.it = TRUE, ...)
{  
  if(is.null(col)) col=1

  # extract parameters of CUTT distribution from y as a list 
  ft = getCUTTparams(y)
  
  # construct function to transform from standard normal values to fitted CUTT distribution
  trn = function(z) qcutt(pnorm(z),theta=ft$theta,K=ft$K,sd=ft$sd,df=ft$df)
  
  # map ages across to normal quantiles using fitted CUTT distribution
  ageNorm = qnorm(pcutt(ft$ages,theta=ft$theta,K=ft$K,sd=ft$sd,df=ft$df))
  
  # simulate normal quantiles to compare them to
  resids = matrix(rnorm(n.sim*length(ft$ages)),ncol=n.sim)
  
  # construct qqenvelope. Suppress warnings because resids and ft$sd have different lengths
  out=suppressWarnings(ecostats:::qqnormEnvelope(ageNorm, resids, n.obs=length(ft$ages), transform=trn, conf.level=conf.level, type=type,
                            plot.it=plot.it, main=c(" ", main), ylab=c(" ", ylab), xlab=c(" ", xlab), col=col, line.col=rep(line.col,2), envelope.col=rep(envelope.col,2), ...))
  if(plot.it)
    invisible(out)
  else
    return(out)
}

