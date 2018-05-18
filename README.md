# HMSC

## Installation notes
Currently the recommended way to install the **Hmsc** package is to use `install_github` function from **devtools** package, available from CRAN. Additionally, since `Hmsc` is dependent on the `BayesLogit`package to support observation models for count data, and `BayesLogit` was temporary removed from CRAN, it is necessary to install it before installing `Hmsc`.

The following lines should sucessfully install Hmsc to your R in most cases:

```R
library(devtools)
install_url('https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz')
install_github("gtikhonov/HMSC")
```
