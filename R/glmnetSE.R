#####################
# glmnetSE function #
#####################

#' @title glmnetSE: Add bootstrap SE to glmnet for selected coefficients (no shrinkage)
#' @description glmnet or cv.glmnet with bootstrap standard errors for selected coefficients with no shrinkage applied for them.
#' @param data A data frame, tibble or matrix object with the outcome variable in the first column and the feature variables in the following columns.
#' @param cf.no.shrnkg A character string of the coefficients that will be interpreted, the inference statistic is of interest and therefore no shrinkage will be applied.
#' @param alpha Alpha value [0,1]. An alpha of 0 results in a ridge regression and a value of 1 in a LASSO.
#' @param method A character string defining if no 10-fold cross validation "none", 10-fold cross validation and selecting the lambda at which the MSE is within one standard error of the smallest MSE  "10CVoneSE", or 10-fold cross validation and selecting the lambda at which the smallest MSE is achieved "10CVmin" - default is "none.
#' @param r Number of bootstrap repetitions -default is 100.
#' @param nlambda Number of tested lambda values - default is 100.
#' @param seed Seed set for the bootstrap sampling - default 0 which means no seed set.
#' @param family A character string representing the used model family either "gaussian" or "binomial" - default is "gaussian".
#' @param type A character string indicating the type of calculated bootstrap intervals. It can be "norm", "basic", "stud", "perc", "bca". For more information check the "boot.ci" package - default is "basic".
#' @param conf Indicates the confidence interval level - default is 0.95.
#' @param perf.metric A character string indicating the used performance metric to evaluate the performance of different lambdas and the final model. Can be either "mse" (mean squared error) or "mae" (mean absolute error). Is not applied when method "none" is used - default is "mse".
#' @return A data frame with all coefficients, confidence intervals and information about significance of the effect for the coefficients of interesst.
#' @keywords glmnet standard errors bootstrap shrinkage
#' @examples
#' \dontrun{
#' # GAUSSIAN with no cross validation and coefficient of interest is Education
#' glmnetSE(data=swiss,cf.no.shrnkg = c("Education"), alpha=1, method="none", r=100, seed = 123, family="gaussian")
#'
#'
#' # BINOMIAL with 10-fold cross validation selecting the lambda at which the
#' # smallest MSE is achieved, 500 bootstrap repetitions and coefficient of interest
#' # are Education and Catholic.
#' glmnetSE(data=swiss,cf.no.shrnkg = c("Education", "Catholic"), alpha=1, method="10CVmin", r=500, seed = 123, family="binomial")
#' }
#' @export

glmnetSE <- function(data, cf.no.shrnkg, alpha=1, method="none", r=100, nlambda=100, seed=0, family="gaussian", type="basic", conf=0.95, perf.metric="mse"){

  # Needed packages
  packages = c("boot", "glmnet")

  ## Now load or install packages
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )


  coef.name <- colnames(data)[-1]
  c <- which(names(data) %in% coef.name)-1
  c <- 1:(ncol(data)-1)
  no.shrink <- which(names(data) %in% cf.no.shrnkg)-1

  p.fac <- rep(1, ncol(data)-1)
  p.fac[no.shrink] <- 0


  "%ni%" <- Negate("%in%")


  if(method=="10CVoneSE"){
    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = cv.glmnet(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)
      return(coef(fit, s = "lambda.1se")[-1])
    }

    if(seed==0){
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }else{
      set.seed(seed)
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }


  }else if(method=="10CVmin"){
    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = cv.glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)
      coef <- coef(fit, s = "lambda.min")[-1]
      metric <- fit$cvm[fit$lambda == fit$lambda.min]
      return(c(coef, metric))
    }

    if(seed==0){
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }else{
      set.seed(seed)
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }


  }else{
    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac)
      return(coef(fit, s = 0.1)[-1])
    }

    if(seed==0){
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }else{
      set.seed(seed)
      results = boot(data=data, statistic=boot.glmnet, R=r)
    }

  }


  name <- NULL
  coef <- NULL
  SE <- NULL
  KI_low <- NULL
  KI_up <- NULL
  sig <- NULL

  if(method == "none"){
    iter.loop <- 0
    for(C in c){
      iter.loop <- iter.loop + 1
      coef.nm <- coef.name[iter.loop]
      KI <- boot.ci(results, type=type, conf=conf, index=C)[[type]][4:5]
      if(coef.nm %ni% cf.no.shrnkg){
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- NA
        SE <- c(SE, SE_new)
        KI_low_new <- NA
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- NA
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 0
        sig <- c(sig, sig_new)
      }else if((KI[1] > 0 & KI[2] > 0) | (KI[1] < 0 & KI[2] < 0)){
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,sd)[C]
        SE <- c(SE, SE_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 1
        sig <- c(sig, sig_new)
      }else{
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,sd)[C]
        SE <- c(SE, SE_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 0
        sig <- c(sig, sig_new)
      }
      output <- list(coefficients = name,
                     variables_no_shrinkage = cf.no.shrnkg,
                     estimat = coef,
                     standard_error = SE,
                     CI_low = KI_low,
                     CI_up = KI_up,
                     significance = sig,
                     metric_name = "not applied")
    }

  }else{
    metric.KI <- boot.ci(results, type=type, conf=conf, index=max(c)+1)[[type]][4:5]
    metric <- results[["t0"]][max(c)+1]
    iter.loop <- 0
    for(C in c){
      iter.loop <- iter.loop + 1
      coef.nm <- coef.name[iter.loop]
      KI <- boot.ci(results, type=type, conf=conf, index=C)[[type]][4:5]
      if(coef.nm %ni% cf.no.shrnkg){
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- NA
        SE <- c(SE, SE_new)
        KI_low_new <- NA
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- NA
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 0
        sig <- c(sig, sig_new)
      }else if((KI[1] > 0 & KI[2] > 0) | (KI[1] < 0 & KI[2] < 0)){
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,sd)[C]
        SE <- c(SE, SE_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 1
        sig <- c(sig, sig_new)
      }else{
        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,sd)[C]
        SE <- c(SE, SE_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        sig_new <- 0
        sig <- c(sig, sig_new)
      }
      output <- list(coefficients = name,
                     variables_no_shrinkage = cf.no.shrnkg,
                     estimat = coef,
                     standard_error = SE,
                     CI_low = KI_low,
                     CI_up = KI_up,
                     significance = sig,
                     metric_name = perf.metric,
                     metric = metric,
                     metric_CI = metric.KI)
    }

  }

  class(output) <- c("glmnetSE")

  invisible(output)

}



#####################
# summary function #
#####################

#' @title Summary Function for glmnetSE Objects
#' @description Print the coefficients with standard errors and confidence intervals of a glmnetSE object. The inference statistics are only available for the coefficients without shrinkage, because otherwise they are biased. If cross fold validation is applied the performance metric can be displayed.
#' @param object A object of the class glmnetSE.
#' @param output A character string indicating if the coefficients with inference statistics "coef", or the performance metric with confidence intervalls should be displayed "metric" - default is "coef".
#' @return The output of a glmnetSE object or the performance metric if cross fold validation is used.
#' @keywords glmnetSE summary results output
#' @examples
#' \dontrun{
#' # Estimate modell
#' glmnetSE.model <- glmnetSE(data=swiss,cf.no.shrnkg = c("Education"))
#'
#' # OUTPUT coefficients
#' summary.glmnetSE(glmnetSE.model, output="coef")
#'
#' # OUTPUT performance metric
#' summary.glmnetSE(glmnetSE.model, output="metric")
#' }
#' @export

summary.glmnetSE <- function(object, output="coef"){

  if(object$metric_name == "not applied"){
    coef <- data.frame(cbind("Coefficients" = object$coefficients,
                             "Estimates" = round(object$estimat,4),
                             "Std.Error" = round(object$standard_error,4),
                             "CI.low" = round(object$CI_low,4),
                             "CI.up" = round(object$CI_up,4),
                             "Sig." = ifelse(object$significance==1, "*", "")))

    metric <- warning("Metric only applicable with method '10CVmin' or '10CVoneSE'")

  }else{
    coef <- data.frame(cbind("Coefficients" = object$coefficients,
                             "Estimates" = round(object$estimat,4),
                             "Std.Error" = round(object$standard_error,4),
                             "CI.low" = round(object$CI_low,4),
                             "CI.up" = round(object$CI_up,4),
                             "Sig." = ifelse(object$significance==1, "*", "")))

    metric <- data.frame(cbind("Metric" = object$metric_name,
                               "Value" = round(object$metric,2),
                               "Intervall.low" = round(object$metric_CI[1],2),
                               "Intervall.up" = round(object$metric_CI[2],2)))
  }

  if(output == "coef"){
    output.df <- coef
  }else{
    output.df <- metric
  }

 return(output.df)
}

