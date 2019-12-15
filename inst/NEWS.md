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

### New Features

* `Hmsc` models can be modified using `update(<Hmsc model>, <new
  arguments>)`. This was achieved by adding a call component per wish
  in [github issue #34](https://github.com/hmsc-r/HMSC/issues/34).

* `evaluateModelFit` can handle probit models where binary data
  were given as `TRUE`/`FALSE`. Earlier only numeric data (`0`/`1`)
  were accepted. See
  [github issue #30](https://github.com/hmsc-r/HMSC/issues/30).

### Documentation

* Added section on priors in vignette on high-dimensional multivariate models.

Version 3.0-2
=============

* This is the first CRAN release. For previous development, see
  https://github.com/hmsc-r/HMSC/.
