# HMSC

## Background

## Installation notes
Currently the recommended way to install the `HMSC` package is to use `install_github` function from the `devtools` package, available from CRAN. Additionally, since `HMSC` is dependent on the `BayesLogit`package to support observation models for count data, and `BayesLogit` was temporary removed from CRAN, it is necessary to install it before installing `HMSC`.

The following lines should sucessfully install Hmsc to your R in most cases:

```R
library(devtools)
install_url('https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz')
install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual")))
```

## Literature
During the development of `Hmsc` several papers have been publisched describing the different components of the model. 
