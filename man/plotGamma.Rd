% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGamma.R
\name{plotGamma}
\alias{plotGamma}
\title{plotGamma}
\usage{
plotGamma(
  hM,
  post,
  param = "Support",
  trOrder = "Original",
  trVector = NULL,
  covOrder = "Original",
  covVector = NULL,
  trNamesNumbers = c(TRUE, TRUE),
  covNamesNumbers = c(TRUE, TRUE),
  supportLevel = 0.9,
  main = NULL,
  cex = c(0.8, 0.8, 0.8),
  colors = colorRampPalette(c("blue", "white", "red")),
  colorLevels = NULL,
  mar = c(6, 9, 2, 0),
  smallplot = NULL,
  bigplot = NULL,
  newplot = TRUE
)
}
\arguments{
\item{hM}{a fitted \code{Hmsc} model object}

\item{post}{posterior summary of Gamma parameters obtained from \code{\link{getPostEstimate}}}

\item{param}{controls which parameter is plotted, current options include "Mean" for posterior mean
estimate, "Support" for the level of statistical support measured by posterior probability for a
positive or negative response, and "Sign" to indicate whether the response is positive,
negative, or neither of these given the chosen \code{supportLevel}}

\item{trOrder}{controls the ordering of traits, current options are "Original", and "Vector".
If trOrder = "Vector", an ordering vector must be provided (see trVector)}

\item{trVector}{controls the ordering of traits if trOrder = "Vector". If a subset of traits
are listed, only those will be plotted}

\item{covOrder}{controls the ordering of covariates, current options are "Original" and
"Vector". If covOrder = "Vector", an ordering vector must be provided (see covVector)}

\item{covVector}{controls the ordering of covariates if covOrder = "Vector". If a subset of
covariates are listed, only those will be plotted}

\item{trNamesNumbers}{logical of length 2, where first entry controls whether trait names
are added to axes, and second entry controls whether traits numbers are added}

\item{covNamesNumbers}{logical of length 2, where first entry controls whether covariate names
are added to axes, and second entry controls whether covariate numbers are added}

\item{supportLevel}{controls threshold posterior support for plotting}

\item{main}{main title for the plot}

\item{cex}{controls character expansion (font size). Three values, controlling covariate names,
trait names, and color legend axis labels}

\item{colors}{controls the colors of the heatmap, default value \code{colorRampPalette(c("blue","white","red"))}}

\item{colorLevels}{number of color levels used in the heatmap}

\item{mar}{plotting margins}

\item{smallplot}{passed to \code{\link[fields]{image.plot}}}

\item{bigplot}{passed to \code{\link[fields]{image.plot}}}

\item{newplot}{set to  false if the plot will be part of multi-panel plot}
}
\description{
Plots heatmaps of parameter estimates or posterior support values of trait effects
on species' environmental responses, i.e. how environmental responses in \code{Beta} responds to
covariates in \code{X}
}
\examples{
# Plot posterior support values of trait effects on environmental responses
gammaPost=getPostEstimate(TD$m, "Gamma")
plotGamma(TD$m, post=gammaPost, param="Support")

# Plot parameter estimates of trait effects on environmental responses
gammaPost=getPostEstimate(TD$m, "Gamma")
plotGamma(TD$m, post=gammaPost, param="Mean")

}
