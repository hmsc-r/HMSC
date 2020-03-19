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
