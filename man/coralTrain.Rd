% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coralTrain.R
\name{coralTrain}
\alias{coralTrain}
\title{coralTrain}
\usage{
coralTrain(Y, XData, XFormula, prior.mu, prior.V, transient, samples, thin)
}
\arguments{
\item{Y}{covariates data.frame}

\item{XData}{arg2}

\item{XFormula}{arg3}

\item{prior.mu}{arg4}

\item{prior.V}{arg5}

\item{transient}{arg6}

\item{samples}{arg7}

\item{thin}{arg8}
}
\value{
list with means and covariance matrices of CORAL posteriors
}
\description{
Trains CORAL models
}
