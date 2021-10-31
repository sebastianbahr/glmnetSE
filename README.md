
# glmnetSE

<!-- badges: start -->
[![R-CMD-check](https://github.com/sebastianbahr/glmnetSE/workflows/R-CMD-check/badge.svg)](https://github.com/sebastianbahr/glmnetSE/actions)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
<!-- badges: end -->

The package glmnetSE allows the user to build a LASSO, Ridge, or Elastic
Net model with bootstrap inference statistics (SE, CI, and p-value) for
selected coefficients with no shrinkage applied for them. The models are
fitted with glmnet or cv.glmnet if cross-validation is desired. If a
test data set is supplied the model performance can be evaluated on the
train as the test set. If an Elastic Net should be fitted it is possible
to pass a sequence of values for alpha. The alpha which results in the
best performance metric is selected for the model fit. Because bootstrap
sampling is time-consuming clusters are built for parallelized
computation. The build-in function summary() allows the user to print
the output of the fitted glmnetSE object. The build-in function plot()
allows the user to display a ROC curve if the performance metric AUC is
selected.

## Installation

You can install the released version of glmnetSE from CRAN

``` r
install.packages("glmnetSE")
```

You can install the latest development version from github with:

``` r
# install.packages("devtools")
devtools::install_github("sebastianbahr/glmnetSE")
```

## Example 1

The following example examines if the percentage of population beyond
primary education level affects the fertility rate on municipal level. A
LASSO model (alpha: 1) is fitted and 10-fold cross-validation is used.
The lambda of the least complex model is selected which is within one
standard deviation of the best performing model. The 95% confidence
intervals are calculated using 250 bootstrap repetitions. Because the
outcome variable fertility is continuous, the family gaussian is
selected and the mean squared error is used as performance metric. The
model output is printed by using the summary() function.

``` r
library(glmnetSE)
data("swiss", package = "glmnetSE")

glmnetSE.model <- glmnetSE(data=swiss, cf.no.shrnkg = c("Education"), alpha=1, method="10CVoneSE", r=250, seed = 123, family="gaussian", perf.metric="mse")

summary(glmnetSE.model)
```

    ## LASSO model:
    ## Outcome variable: Fertility 
    ## Variables without shrinkage: Education 
    ## 
    ## 
    ##      Coefficients Estimates Std.Error  CI.low   CI.up p.value Sig.
    ##       Agriculture         0      <NA>    <NA>    <NA>    <NA>     
    ##       Examination         0      <NA>    <NA>    <NA>    <NA>     
    ##         Education   -0.8087    0.1966 -1.2965 -0.4139   0.008   **
    ##          Catholic    0.0507      <NA>    <NA>    <NA>    <NA>     
    ##  Infant.Mortality    0.6573      <NA>    <NA>    <NA>    <NA>     
    ## 
    ## 
    ## Performance Metric: 
    ##  Metric Value CI.low  CI.up
    ##     mse 73.15  58.47 111.87

## Example 2

The following example examines if the percentage of population beyond
primary education level affects the fertility rate on municipal level. A
train and test data set is build to evaluate the model performance on
the test set. A Elastic Net model with a sequence of alphas is fitted
and 10-fold cross-validation is used. The lambda of the least complex
model is selected which is within one standard deviation of the best
performing model. The 95% confidence intervals are calculated using 250
bootstrap repetitions. The outcome variable fertility is dichotomized at
the median. Because the outcome variable is dichotomous, the family
binomial is selected and the AUC is used as performance metric. The
model output is printed by using the summary() function and the ROC
curve on the test data is displayed with the plot() function.

``` r
library(glmnetSE)
data("swiss", package = "glmnetSE")

swiss$Fertility <- ifelse(swiss$Fertility >= median(swiss$Fertility), 1, 0)

set.seed(123)
train_sample <- sample(nrow(swiss), 0.7*nrow(swiss))
swiss.train <- swiss[train_sample, ]
swiss.test  <- swiss[-train_sample, ]

glmnetSE.model <- glmnetSE(data=swiss.train, cf.no.shrnkg = c("Education"), alpha=seq(0.1,0.9,0.1), method="10CVoneSE", test=swiss.test, r=250, seed = 123, family="binomial", perf.metric="auc")
```

    ## [1] "The best alpha is:  0.1"

``` r
summary(glmnetSE.model)
```

    ## Elastic Net model:
    ## Outcome variable: Fertility 
    ## Variables without shrinkage: Education 
    ## 
    ## 
    ##      Coefficients Estimates Std.Error  CI.low  CI.up p.value Sig.
    ##       Agriculture   -0.0032      <NA>    <NA>   <NA>    <NA>     
    ##       Examination   -0.0308      <NA>    <NA>   <NA>    <NA>     
    ##         Education   -0.1692    0.2382 -0.2943 0.7264     0.2     
    ##          Catholic    0.0099      <NA>    <NA>   <NA>    <NA>     
    ##  Infant.Mortality      0.09      <NA>    <NA>   <NA>    <NA>     
    ## 
    ## 
    ## Performance metric train data: 
    ##  Metric Value CI.low CI.up
    ##     auc   0.9   0.84     1
    ## 
    ## 
    ## Performance metric test data: 
    ##  Metric Value CI.low CI.up
    ##     auc   0.7   0.43     1

``` r
plot(glmnetSE.model)
```

![](README_files/figure-gfm/example2-1.png)<!-- -->
