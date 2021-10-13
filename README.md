
# glmnetSE

<!-- badges: start -->
[![R-CMD-check](https://github.com/sebastianbahr/glmnetSE/workflows/R-CMD-check/badge.svg)](https://github.com/sebastianbahr/glmnetSE/actions)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<!-- badges: end -->

The package "glmnetSE" allows the user to obtain bootstrap standard errors and confidence intervals for LASSO, ridge and elastic net models of the "glmnet" package. Because these regularization methods are implementing a shrinkage to the coefficients they represent not the true effect size and should neither be interpreted nor the estimation of standards errors is advisable. The package "glmnetSE" applies no shrinkage parameter to selected variables of interest, which than can be interpreted by the user. For the selected variables bootstrap standard errors and confidence intervals are estimated. The implementation of "cv.glmnet" gives the user the ability to build robust models with 10-fold cross validation. Further, it is possible to select the least complex model within one standard error of the smallest model error and the model with the smallest model error. As model error MSE and MAE can be used.   

## Installation

You can install the released version of glmnetSE from [CRAN]
``` r
install.packages("glmnetSE")
```

You can install the latest development version from github with:
```r
# install.packages("devtools")
devtools::install_github("sebastianbahr/glmnetSE")
```


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(glmnetSE)
glmnetSE.model <- glmnetSE(data=swiss,cf.no.shrnkg = c("Education"), alpha=1, method="none", r=100, seed = 123, family="gaussian")

summary(glmnetSE.model)
```

