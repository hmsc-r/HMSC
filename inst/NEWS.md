Version 3.0-13
==============

### New Features

* Posterior samples from several independent `Hmsc` objects can be
  combined as new chains with new method function `c()`. This provides
  an easy alternative for distributed computing. The user should take
  care that these independent models are defined similarly so that
  they really can be combined. The function tests for similarity of
  objects, but it only gives warnings and can allow combination of
  incompatible models at user will. The user should be careful **not**
  to start these models from the same random number seed as these just
  duplicate your data instead of adding independent new samples.

* `sampleMcmc` allows the use of fork clusters instead of socket
  clusters. Socket clusters are still the only alternative in Windows,
  but other platforms can profit from the use of fork clusters which
  may have lower memory use and are faster to set up, and also may be
  marginally faster. The choice can be made with new argument
  `useSocket` (defaults `TRUE`).

* Updaters in `sampleMcmc` can occasionally fail in extreme `Hmsc`
  models. This is no longer an error that stops analysis, but sampling
  tries to recover from failures. The numbers of failures for each
  updater is reported with the result. If there are only a low number
  of failures, the sampling is safe to use. If there are updater
  errors only in some chains, these chains can be removed, but other
  chains can be used. See
  [issue #123](https://github.com/hmsc-r/HMSC/issues/123).

* New experimental function `pcomputePredictedValues` with more
  aggressive parallelization than `computePredictedValues`. In old
  code chains within each partition could be run in parallel, but
  partitions were run serially. In the new function, all chains and
  partitions can be run in parallel. The plan is to replace the old
  function with this new alternative, but at the moment both functions
  are available for testing. See
  [issue #142](https://github.com/hmsc-r/HMSC/issues/142).

* Implemented longitude-latitude coordinates and user-supplied
  distance matrices for NNGP spatial models. Sanity checks for spatial
  model input were improved.

* Improved support for spatial models defined _via_ distance matrices
  instead of spatial coordinates.

* `constructGradient` provides wider choice of coordinates for
  centroid of `new_unit`, including user-set and infinite (meaning no
  spatial dependence) coordinates.

* Detect cases when user tries to analyse posterior samples of
  non-sampled `Hmsc` object to avoid confusing error messages such as
  reported in [issue #125](https://github.com/hmsc-r/HMSC/issues/125).

### Bug Fixes

* `sampleMcmc` with `initPar = "fixed effects"` failed if **Y**
  variates had missing values. The choice `"fixed effects"` was
  undocumented in the package, but was used in several scripts at
  large. See [issue #101](https://github.com/hmsc-r/HMSC/issues/101).

* The default number of neighbours in NNGP spatial models was not
  known in all posterior analysis tools giving very obscure error
  messages. Reported by Ben Weigel (Uni Helsinki).

* Covariate-dependent latent loadings did not have correct alignment.

* `predict` did not honour setting `start` and `thin` which could
  result in huge output data that exhausted memory. See
  [issue #86](https://github.com/hmsc-r/HMSC/issues/86).

* `predict` failed in NNGP spatial models. See
  [issue #96](https://github.com/hmsc-r/HMSC/issues/96).

* `predict` failed with one-dimensional spatial data. See comments in
  unrelated [issue #61](https://github.com/hmsc-r/HMSC/issues/61).

* Missing values are handled better in `predict`, but they are still
  not allowed in all cases.

* `prepareGradient` failed with geo-referenced spatial random levels.

* More robust handling of models that were fitted with model matrix
  `X` instead of model frame `XData` and model formula
  `XFormula`. Concerns functions `biPlot`,
  `computeVariancePartitioning` and `constructGradient`. Fixes
  [issue #126](https://github.com/hmsc-r/HMSC/issues/126).

* `biPlot` has improved handling of colour scaling of continuous
  variables.

Version 3.0-11
==============

### Installation

* New package dependencies: **matrixStats** and **pracma**.

### New Features

* `plotGradient` gained argument to show the support of trend for
  continuous variables. Main title can be shown for factor variables
  (earlier it was shown only for continuous variables).

* Grids of knots for Gaussian Predictive Process (GPP) are centred for
  the coordinates in `constructKnots`. More knots were produced than
  requested.

### Bug Fixes

* Prediction failed in spatial models with `predictEtaMean = TRUE`.

* Prediction failed in spatial NNGP models.

* `constructGradient` (and hence `plotGradient`) ignored specified
  order of factor levels. See github
  [issue #63](https://github.com/hmsc-r/HMSC/issues/63).

* Performance inefficiency issues were fixed in NNGP models and some
  updaters.

### User Interface

* User interface is more robust and can handle several inputs that
  earlier caused errors (often with confusing and obscure error
  messages). Input data is checked more carefully to avoid misleading
  results because of wrongly interpreted data.

* User interface changes fix several github issues:
  [#65](https://github.com/hmsc-r/HMSC/issues/65),
  [#66](https://github.com/hmsc-r/HMSC/issues/66),
  [#68](https://github.com/hmsc-r/HMSC/issues/68),
  [#70](https://github.com/hmsc-r/HMSC/issues/70),
  [#71](https://github.com/hmsc-r/HMSC/issues/71),
  [#78](https://github.com/hmsc-r/HMSC/issues/78),
  [#80](https://github.com/hmsc-r/HMSC/issues/80),
  [#81](https://github.com/hmsc-r/HMSC/issues/81),
  [#82](https://github.com/hmsc-r/HMSC/issues/82).

* Spatial and phylogenetic data are inspected more carefully to avoid
  errors in sampling.

* Updaters are automatically disabled when needed instead of producing
  an error.
 
Version 3.0-9
=============

### Installation

* **Hmsc** is no longer dependent on packages **mvtnorm** and
  **pdist**.

* **Hmsc** is now dependent on the **sp** package.

* Vignettes can be re-built from their sources out of the box.
  Previously they needed editing by hand to reproduce the pdf version.

### Spatial Data

* It is now possible to use Spatial data in random models. Handling of
  Spatial data is based on the **sp** package and follows its
  conventions. The locations of sampling units can be given as a
  decimal longitude-latitude matrix, and the **Hmsc** functions will
  use great circle distances in spatial models. Projected spatial
  coordinates will be handled as such and Euclidean distances will be
  used internally.

* User-specified spatial distances can be more widely used in spatial
  random models. However, some models are more flexible with spatial
  coordinates. Most importantly, Gaussian Predictive Process (GPP)
  needs spatial coordinate data.

### New Features

* Species data `Y` is normally a numeric matrix, but now it is allowed
  to use numeric data frames, or in univariate models, a numeric vector.

* A `tibble` can be used for measured covariates for fixed effects
  `XData` in addition to a data frame (the wish of Github
  [issue #37](https://github.com/hmsc-r/HMSC/issues/37)).

* The names of `distr`ibutions can be abbreviated in `Hmsc` definition
  as long as the names are unique.

* `computeWAIC` is more robust against results of poorly fitting
  models, and it is now possible to evaluate WAIC separately for each
  species. See GitHub
  [issue #44](https://github.com/hmsc-r/HMSC/issues/44).

* `constructGradient` argument `nonFocalVariables` accepts now a
  single number `1` or `2` as a shortcut of default type for all
  non-focal variables instead of requesting a list of types of all
  variables.

* `plotGradient` gained new argument `yshow` which is a single number
  or vector of numeric values that must be included on the
  *y*-axis. In general, the *y*-axis is scaled to show the plotted
  values, but `yshow = 0` will always show zero, even when this is not
  among plotted values, and `yshow = c(0,1)` will show both zero and
  one.

* `plotVariancePartition` defaults to plot the original terms instead
  of single contrast. For instance, only one component is shown for
  multilevel factors instead of showing each level separately. User
  can still specify how the components are displayed.

* plot functions `plotBeta`, `plotGamma` and
  `plotVariancePartitioning` allow setting or modifying the plot main
  title. `plotGradient` already allowed this.

* Random seed is now saved in `sampleMcmc` models. This allows
  replication of same random number sequences. However, there is no
  guarantee of replication across **Hmsc** release versions or
  computing platforms.

* `HmscRandomLevel` saves the function call. The call can be inspected
  with `getCall()` and the model can be modified with `update()`.


### Bug Fixes

* `constructGradient` could sometimes shuffle spatial locations
  leading into wrong predictions with spatial models.

* `plotGradient(..., showData = TRUE)` ignored data values in setting
  plot minimum. See GitHub
  [issue #48](https://github.com/hmsc-r/HMSC/issues/48). The data
  values were not always shown with `measure = "S"` in quantitative
  linear models.

Version 3.0-6
=============

### Installation and Versions

* **R** release 4.0 will drop the convention to automatically change
  character variables to factors, and this causes errors in internal
  working of several **Hmsc** functions. This version of **Hmsc** is
  released principally to accomodate these changes in **R**. **Hmsc**
  will also work in previous versions of **R**.

* **Hmsc** 3.0-5 was never released to CRAN. It is a snapshot that
  corresponds to the on-line publication of Tikhonov *et al.* (2020)
  Joint species distribution modelling with the R-package
  Hmsc. *Methods in Ecology and Evolution* **11,**
  442--447. (https://doi.org/10.1111/2041-210X.13345).

### Models

* Shape and rate parameters (`aSigma`, `bSigma`) for the prior Gamma
  distribution for the variance parameter (`sigma`) changed. The
  change will influence models with `"normal"` and `"lognormal
  poisson"` distributions. In particular, `"lognormal poisson"` will
  more easily tend toward zero `sigma` if there is no overdispersion
  to `"poisson"`.  However, in such cases it may be wiser to refit
  models with pure `"poisson"` distribution. You can changes these
  parameters with `setPriors` function.

* Cross-validation works also when the test data set has some spatial
  units that were unseen in the training data.

### Bug Fixes

* When calling `sampleMcmc` with `fromPrior = TRUE`, the residual
  variance parameter `sigma` used Gamma rather than inverse of Gamma
  distribution. The same error was present when sampling the initial
  values for the MCMC algorithm. However, the actual MCMC algorithm
  (and thus the posterior distribution) was correct.

* Predictions with spatial NNGP models failed if there was only one
  unit. Github [issue #40](https://github.com/hmsc-r/HMSC/issues/40).

* Reduced-Rank Regression also works for single-species models, and
  more robust scaling is used for species-specific covariate matrices.
  
* Spatial models with Gaussian Predictive Process now also works when
  the number of spatial locations is less than the number of sampling 
  units.

* Predictions with spatial NNGP and GPP models gave bad estimates. 

### Documentation

* New vignette (#5) on **Hmsc** performance.

Version 3.0-4
=============

### Installation

* **Hmsc** no longer depends on **phytools** package, and external
  **ImageMagick** software is no longer needed. See discussion in
  [github issue #34](https://github.com/hmsc-r/HMSC/issues/34) and
  [pull request #36](https://github.com/hmsc-r/HMSC/pull/36).

### Bug Fixes

* Several functions failed in the development version of **R** (to be
  released as **R** version 4). The failures were caused by changes in
  **R** internals.

* Fixed bug with delta for `alignPosterior` which influences
`sampleMcmc`. See
[github issue #27](https://github.com/hmsc-r/HMSC/issues/27).

* `plotBeta` failed with argument `plotTree = FALSE` together with
`SpeciesOrder = "Tree"`.

* Spatial models with Nearest Neighbour Gaussian Process (NNGP) failed
  when the number spatial locations was not equal to the number of
  sampling units. This could happen, for instance, if there are
  multiple observations on the same spatial location. The problem
  still persists in spatial models with Gaussian Predictive Process
  (GPP).

### New Features

* `Hmsc` models can be modified using `update(<Hmsc model>, <new
  arguments>)`. This was achieved by adding a call component per wish
  in [github issue #34](https://github.com/hmsc-r/HMSC/issues/34).

* `evaluateModelFit` can handle probit models where binary data
  were given as `TRUE`/`FALSE`. Earlier only numeric data (`0`/`1`)
  were accepted. See
  [github issue #30](https://github.com/hmsc-r/HMSC/issues/30).

* `biPlot` uses equal aspect ratio in ordination biplots.

### Documentation

* Added section on priors in vignette on high-dimensional multivariate models.

Version 3.0-2
=============

* This is the first CRAN release. For previous development, see
  https://github.com/hmsc-r/HMSC/.
