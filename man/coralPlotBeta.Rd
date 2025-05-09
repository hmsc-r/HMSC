% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coralPlotBeta.R
\name{coralPlotBeta}
\alias{coralPlotBeta}
\title{coralPlotBeta}
\usage{
coralPlotBeta(
  mu,
  V,
  phyloTree,
  spNames.common,
  plotColumns = c(1:ncol(mu)),
  quantile.support = c(0.05, 0.95),
  plotTree = 0.3,
  seed = NULL,
  col.common = c("blue", "red"),
  alpha.common = 0.5,
  jitter.common = 0.45,
  cex.common = 0.5,
  pch.common = 16,
  col.rare = c("cyan", "pink"),
  alpha.rare = 0.5,
  jitter.rare = 0.2,
  cex.rare = 0.2,
  pch.rare = 16,
  showCovNames = TRUE
)
}
\arguments{
\item{m}{fitted \code{Hmsc}-class object}

\item{muList.coral}{arg2}

\item{VList.coral}{arg3}
}
\value{
list with combined means and covariance matrices
}
\description{
Plots summary of Beta coefficents fitted with CORAL approach
}
