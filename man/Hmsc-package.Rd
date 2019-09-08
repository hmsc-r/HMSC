% Maintained manually: no R/Hmsc-package.R file
\docType{package}
\name{Hmsc-package}
\alias{Hmsc-package}
\title{Hmsc: Hierarchical Model of Species Communities}
\description{
Hierarchical Modelling of Species Communities (Hmsc) is a
 flexible framework for Joint Species Distribution Modelling (JSDMs). The
 framework can be used to relate species occurrences or abundances to
 environmental covariates, species traits and phylogenetic relationships.
 JSDMs are a special case of species distribution models (SDMs) that take
 into account the multivariate nature of communities which allows us to
 estimate community level responses as well capture biotic interactions and
 the influence of missing covariates in residual species associations. The
 Hmsc package contains functions to fit JSDMs, analyze the output and to
 generate predictions with these JSDMs.
}
\section{Workflow}{

 A typical workflow for a Hmsc analysis constists of 5 steps:

 \bold{Step 1: Setting the model structure and fitting the model}

 The obligatory data for a Hmsc analysis includes a matrix of species occurrences
 or abundances (Y) and a data frame of environmental covariates (XData). The
 species matrix Y consists of ns columns representing the species and np rows
 representing the sampling units. Species data for different species does not
 have to be on the same scale, i.e. species 1 may be recorded as
 presence/absence while species 2 is recorded as abundance. XData consists of
 nc columns representing the environmental variables and np rows representing
 the sampling units. Y and XData need to have the same amount of rows
 (sampling units).

 Optionally, the user can include information species traits, phylogenetic
 relationships and information on the spatiotemporal or hierarchical context of the sampling
 design to account for dependencies among the sampling units.

 The model structure is created using the \code{\link{Hmsc}} function. As
 input, this function needs at least the matrix with species data (Y) and the
 dataframe with environmental variables. This is also the place where the data model is specified, as a default Hmsc assumes normally distributed species data. Other options are 'probit' for binary data, 'Poisson' and 'overdispersedPoisson' for count data. Additionally, you can specify the study design, the random effects structure, the species traits to be used,
 species phylogeny, how the covariate should be scales and how viable selection should be applied.

 The random levels supplied to \code{\link{Hmsc}} are generated using \code{\link{HmscRandomLevel}}. Here, the user can specify the spatial or temporal information for the units as well as the covariates for covariate-depedent associations. If there is no spatial, temporal or coviate information for this level, the user should provide here the names of the units for this level.

 After setting the model structure, the model is fitted by running \code{\link{sampleMcmc}} to sample from the posterior distributions of the model parameters.

 \bold{Step 2: Examining MCMC convergence}

 After fitting the model, the MCMC convergence needs to be evaluated. The
 easiest way to do this is generate a coda object from the fitted model using
 \code{\link{convertToCodaObject}}. Convergence of the chains should be
 assessed qualitatively using trace plots of the model parameters and
 quantitatively using gelman diagnostics \code{\link[coda]{gelman.diag}} and
 by calculating the effective sample size with
 \code{\link[coda]{effectiveSize}}.

 \bold{Step 3: Evaluating model fit}

 If MCMC convergence is satisfacotry, the model fit can be evaluated. Using
 the \code{\link{computePredictedValues}} function. It is also recommended to
 also compute predictive power using cross validation. Cross validation is
 done by supplying the partitioning created using using
 \code{\link{createPartition}} to \code{\link{computePredictedValues}}.

 \bold{Step 4: Exploring parameter estimates}

 When the fitted model has satisfactory convergence and fit, the next step is
 generally to investigate the parameter estimates of the model. The posterior
 distribution of the parameter of choice is extracted from the model object
 using \code{\link{getPostEstimate}}. \code{\link{plotBeta}} and
 \code{\link{plotGamma}} can be used to visuzalize the beta and gamma
 parameter estimates. Additionally, the function \code{\link{biPlot}} can be
 used to construct an ordination plot from the eta and lamda parameters

 Additionally, at this stage of the analysis, the user may want to look at
 how the variance explained by the model is partitioned using the
 \code{\link{computeVariancePartitioning}} and the
 \code{\link{plotVariancePartitioning}} function.

 \bold{Step 5: Making predictions}

 Hmsc comes with a generic \code{\link{predict}} function that can be used to predict species occurences or counts in new units. These predictions can be unconditional or conditional on the occurence of other species. For conditional predictions, the user needs to supply species data for at least some species in the new units.
 Additionally, the package specific tools to make predictions along environmental gradients. These predictions can be vizualized at both the species and the community level (\code{\link{constructGradient}}, \code{\link{plotGradient}}).

}

\section{Citing the Package}{
 Tikhonov, G., Opedal, O., Abrego, N., Lehikoinen, A., Ovaskainen, O., 2019. Joint species distribution modelling with HMSC-R. \emph{bioRxiv}
}

\seealso{
Useful links:
\itemize{
  \item \url{https://www.helsinki.fi/en/researchgroups/metapopulation-research-centre/hmsc}
  \item Report bugs at \url{https://github.com/hmsc-r/HMSC/issues}
}

}
\keyword{internal}
