% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reginv_S3methods.R
\name{qqenvelope}
\alias{qqenvelope}
\title{Quantile-Quantile Plots with Global Simulation Envelopes}
\usage{
qqenvelope(y, n.sim = 199, conf.level = 0.95, ylab = "Sample Quantiles", ...)
}
\arguments{
\item{y}{is a fitted object.}

\item{n.sim}{the number of simulated sets of residuals to be generated, to which
the observed residuals will be compared. The default is 199 datasets.}

\item{conf.level}{the confidence level to use in constructing the envelope.}

\item{ylab}{\code{y} axis label (if a plot is produced).}

\item{...}{further arguments sent through to \code{plot}.}
}
\description{
Produces a QQ plot from fitted object \code{y} with global envelopes constructed
by simulation. Global envelopes are constructed using the \code{GET} package for 
simultaneous control of error rates over the whole plot. For details see \code{\link[ecostats]{qqenvelope}}
}
