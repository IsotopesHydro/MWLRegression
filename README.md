# MWLRegression

R package to perform linear regression of data with measurement uncertainties in both variables.
The package was developed in the framework of regression of water stable isotopes in precipitation, but can be used for any X and Y dataset with uncertainties in both variables.

# How to Install

You can install this through use devtools:

```r
devtools::install_github("IsotopesHydro/MWLRegression")
```


# Overview

To see the full list of exported functions:

```{r}
library("MWLRegression")
ls("package:MWLRegression")
```

A quick overview of some of the functions:

* `YorkRegression`: Performs York linear regression which accounts for uncertainties in both y and y variable, with uncertainties also potentially correlated. The function depend on 'ggplot2' for plotting. 

* `RegEqualSigma`: Performs linear regression on x and y data which as equal measurement uncertainties. The function depend on 'ggplot2' for plotting.