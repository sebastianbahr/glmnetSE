#####################
# glmnetSE function #
#####################

#' @title glmnetSE: Add Nonparametric Bootstrap SE to Glmnet for Selected Coefficients (No Shrinkage)
#' @description Builds a LASSO, Ridge, or Elastic Net model with \code{\link[glmnet:glmnet]{glmnet}} or \code{\link[glmnet:cv.glmnet]{cv.glmnet}} with bootstrap inference statistics (SE, CI, and p-value) for selected coefficients with no shrinkage applied for them. Model performance can be evaluated on test data and an automated alpha selection is implemented for Elastic Net. Parallelized computation is used to speed up the process.
#' @param data A data frame, tibble, or matrix object with the outcome variable in the first column and the feature variables in the following columns. Note: all columns beside the first one are used as feature variables. Feature selection has to be done beforehand.
#' @param cf.no.shrnkg A character string of the coefficients whose effect size will be interpreted, the inference statistic is of interest and therefore no shrinkage will be applied.
#' @param alpha Alpha value [0,1]. An alpha of 0 results in a ridge regression, a value of 1 in a LASSO, and a value between 0 and 1 in an Elastic Net. If a sequence of possible alphas is passed to the \code{alpha} argument the alpha of the best performing model (based on the selected \code{method} and \code{perf.metric}) is selected - default is 1.
#' @param method A character string defining if 10-fold cross-validation is used or not. Possible methods are \code{none}: no cross-validation is applied and the coefficients for lambda = 0.1 are selected. \code{10CVoneSE }:  10-fold cross-validation is applied and the lambda of the least complex model with an MSE within one standard error of the smallest MSE is selected.  \code{10CVmin}: 10-fold cross-validation is applied and the lambda at which the MSE is the smallest is selected - default is \code{10CVoneSE}.
#' @param test A data frame, tibble, or matrix object with the same outcome and feature variables as supplied to \code{data} which includes test-observations not used for the training of the model.
#' @param r Number of nonparametric bootstrap repetitions - default is 250
#' @param nlambda Number of tested lambda values - default is 100.
#' @param seed Seed set for the cross-validation and bootstrap sampling - default 0 which means no seed set.
#' @param family A character string representing the used model family either \code{gaussian} or \code{binomial} - default is \code{gaussian}.
#' @param type A character string indicating the type of calculated bootstrap intervals. It can be \code{norm}, \code{basic},  \code{perc}, or  \code{bca}. For more information check the \code{\link[boot:boot.ci]{boot.ci}} package - default is \code{basic}.
#' @param conf Indicates the confidence interval level - default is 0.95.
#' @param perf.metric A character string indicating the used performance metric to evaluate the performance of different lambdas and the final model. Can be either \code{mse} (mean squared error), \code{mae} (mean absolute error), \code{class} (classification error), or \code{auc} (area under the curve). Is not applied when method \code{none} is used - default is \code{mse}.
#' @return \code{glmnetSE } object which output can be displayed using \code{summary()} or \code{summary.glmnetSE()}. If family \code{binomial} and performance metric \code{auc} is used it is possible to plot the ROC curve with \code{plot()} or \code{plot.glmnetSE()}.
#' @keywords glmnet standard errors bootstrap shrinkage
#' @author  Sebastian Bahr, \email{sebastian.bahr@@unibe.ch}
#' @references  Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. \url{https://www.jstatsoft.org/v33/i01/}.
#' @references Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software, 39(5), 1-13. \url{https://www.jstatsoft.org/v39/i05/}.
#' @references Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall. \url{https://cds.cern.ch/record/526679/files/0412042312_TOC.pdf}
#' @seealso  \code{\link{summary.glmnetSE}} and \code{\link{plot.glmnetSE}} methods.
#' @examples
#' \dontrun{
#' # LASSO model with gaussian function, no cross validation, a seed of 123, and
#' # the coefficient of interest is Education.
#'
#' glmnetSE(data=swiss, cf.no.shrnkg = c("Education"), alpha=1, method="none", seed = 123)
#'
#'
#' # Ridge model with binomial function, 10-fold cross validation selecting the lambda
#' # at which the smallest MSE is achieved, 500 bootstrap repetitions, no seed, the
#' # performance metric AUC, and the coefficient of interest are Education and Catholic.
#'
#' # Generate dichotom variable
#' swiss$Fertility <- ifelse(swiss$Fertility >= median(swiss$Fertility), 1, 0)
#'
#' glmnetSE(data=swiss, cf.no.shrnkg = c("Education", "Catholic"), alpha=0, method="10CVmin", r=500,
#'          seed = 0, family="binomial", perf.metric = "auc")
#'
#'
#' # Elastic Net with gaussian function, automated alpha selection, selection the lambda
#' # within one standard deviation of the best model, test data to obtain the performance
#' # metric on it, a seed of 123, bias-corrected and accelerated confidence intervals, a
#' # level of 0.9, the performance metric MAE, and the coefficient of interest is Education.
#'
#' # Generate a train and test set before running the code!!!
#'
#' glmnetSE(data=swiss.train, cf.no.shrnkg = c("Education"), alpha=seq(0.1,0.9,0.1),
#' method="10CVoneSE", test = swiss.test, seed = 123, family = "gaussian", type = "bca",
#' conf = 0.9, perf.metric = "mae")
#' }
#' @export

glmnetSE <- function(data, cf.no.shrnkg, alpha=1, method="10CVoneSE", test="none",r=250, nlambda=100, seed=0, family="gaussian", type="basic", conf=0.95, perf.metric="mse"){

  # Needed packages
  packages = c("boot", "glmnet", "parallel")

  # Install packages
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        utils::install.packages(x, dependencies = TRUE)
      }
    }
  )

  # Get coefficient names and amount
  coef.name <- colnames(data)[-1]
  c <- which(names(data) %in% coef.name)-1
  c <- 1:(ncol(data)-1)
  n <- nrow(data)
  # Get number coefficients without shrinkage
  no.shrink <- which(names(data) %in% cf.no.shrnkg)-1

  # Build penalty factor
  p.fac <- rep(1, ncol(data)-1)
  p.fac[no.shrink] <- 0

  # Get name outcome variable
  dep.var <- colnames(data[1])
  options(warn=-1)

  # Get model name
  if(alpha == 0){

    model_type = "Ridge"

  }else if(alpha == 1){

    model_type = "LASSO"

  }else{

    model_type = "Elastic Net"

  }

  options(warn = 1)
  "%ni%" <- Negate("%in%")

  # Warnings section

  # Warning and stop if data input is not a data frame, tibble or matrix
  if(class(data)!= "data.frame" && class(data)!= "tbl_df" && class(data)!= "matrix"){
    warning("Use an object of class data frame, tibble or matrix as data input")
    stop()
  }

  # Warning and stop if cf.no.shrnkg are not character values
  if(class(class(cf.no.shrnkg))!= "character"){
    warning("The input of the coefficients without applied shrinkage has to be a character")
    stop()
  }

  # Warning and stop if alpha is not numeric and not between 0 and 1
  if(class(alpha)!= "numeric"){
    warning("Alpha has to be numeric")
    stop()
  }else if(any(alpha > 1) || any(alpha < 0)){
    warning("Alpha should be between 0 and 1")
    stop()
  }

  # Warning and stop if method is not 'none', '10CVmin' or '10CVoneSE'
  if(method!= "none" && method!= "10CVmin" && method!= "10CVoneSE"){
    warning("The method should be 'none', '10CVmin' or '10CVoneSE'")
    stop()
  }

  # Warning and stop if bootstrap repetitions are not numeric and just warning if they smaller than 100.
  if(class(r)!= "numeric"){
    warning("The number of bootstrap repetitions has to be numeric")
    stop()
  }else if(r < 100){
    warning("Be aware that you are using less than 100 bootstrap repetitions")
  }

  # Warning and stop if tested lambdas are not numeric.
  if(class(nlambda)!= "numeric"){
    warning("The number of tested lambdas has to be numeric")
    stop()
  }

  # Warning and stop if model family is not 'gaussian' or 'binomial'
  if(family!= "gaussian" && family!= "binomial"){
    warning("The model family should be 'gaussian' or 'binomial'")
    stop()
  }

  # Warning and stop if bootstrap CI type is not 'basic', 'norm', 'perc', or 'bca'.
  if(type!= "basic" && type!= "norm" && type!= "perc" && type!= "bca"){
    warning("The bootstrap confidence interval type has to be 'basic', 'norm', 'perc', or 'bca'")
    stop()
  }

  # Warning and stop if significance level is not numeric and not between 0 and 1
  if(class(conf)!= "numeric"){
    warning("Conf has to be numeric")
    stop()
  }else if(any(conf > 1) || any(conf < 0)){
    warning("Conf should be between 0 and 1")
    stop()
  }

  # Warning and stop if performance metric is not 'mse', 'mae', 'auc' or 'class'
  if(perf.metric!= "mse" && perf.metric!= "mae" && perf.metric!= "auc" && perf.metric!= "class"){
    warning("The performance metric should be 'mse', 'mae', 'class' or 'auc'")
    stop()
  }

  # Warning and stop if family 'binomial' and not performance metric 'auc' or 'class' are used.
  if(family == "binomial" && (perf.metric!= "auc" && perf.metric!= "class")){
    warning("With family 'binomial' please use the performance metric 'class' or 'auc'")
    stop()
  }

  # Warning and stop if family 'gaussian' and not performance metric 'mse' or 'mae' are used.
  if(family == "gaussian" && (perf.metric!= "mse" && perf.metric!= "mae")){
    warning("With family 'gaussian' please use the performance metric 'mse' or 'mae'")
    stop()
  }

  # Warning if test data supplied but not method 10CVmin od 10CVoneSE used
  if((class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix") & (method != "10CVmin" && method != "10CVoneSE")){
    warning("The model performance on the test data is only supplied for method '10CVmin' and '10CVoneSE'")
  }

  # Warning and stop if bootstrap repetitions are smaller than sample size. Ratio important for "bca"
  if((r < n) && type=="bca"){
    warning("The number of bootstrap repetitions 'r' should at leaste be as large as the sample size.")
    stop()
  }

  # Warning if automated alpha selection used without setting a seed
  if((length(alpha)>1) && seed == 0){
    warning("For using automated alpha selection properly please set a seed. The seed should not be 0.")
  }

  # Get best performing alpha value for Elastic Nets
  if(length(alpha)>1){

    if(all(alpha != 0) && all(alpha != 1)){

      enet <- NULL

      for(i in alpha){

        y = as.matrix(data[,1])
        x = as.matrix(data[,-1])
        set.seed(seed)
        enet.fit <- glmnet::cv.glmnet(x=x, y=y, alpha=i, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)
        enet.fit[["alpha"]] = i
        enet.fit_new = enet.fit
        enet[[length(enet) + 1]] = enet.fit_new

      }

      alpha.n <- 1:length(alpha)
      enet.mod.check <- data.frame(alpha = numeric(), mtr.min = numeric(),
                                   mtr.oneSE = numeric())
      for(i in alpha.n){

        mtr.min = enet[[i]]$cvm[enet[[i]]$lambda == enet[[i]]$lambda.min]
        mtr.oneSE = enet[[i]]$cvm[enet[[i]]$lambda == enet[[i]]$lambda.1se]
        alpha = enet[[i]]$alpha
        enet.mod.check[nrow(enet.mod.check) + 1,] = c(alpha, mtr.min, mtr.oneSE)

      }

      if(method == "10CVmin"){

        alpha <- enet.mod.check$alpha[enet.mod.check$mtr.min==min(enet.mod.check$mtr.min)]
        print(paste("The best alpha is: ", alpha))

      }else if(method == "10CVoneSE"){

        alpha <- enet.mod.check$alpha[enet.mod.check$mtr.oneSE==min(enet.mod.check$mtr.oneSE)]
        print(paste("The best alpha is: ", alpha))

      }else{

        warning("Select method '10CVmin' or '10CVoneSE'")

      }
    }else{

      warning("A Elastic Net should have a alpha >0 and <1. Alpha values of 0 and 1 results in a ridge or LASSO model.")

    }
  }else{

    alpha <- alpha

  }


  # Function for AUC and ROC
  getROC_AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
  }


  # Set everything up for cluster computation with maximal use of 32 CPUs because efficiency starts to drop (Sloan et al. 2014)
  envir = environment(glmnetSE)
  cpus <- parallel::detectCores()
  cl <- parallel::makeCluster(ifelse((cpus-1)>32, 32, (cpus-1)))
  parallel::clusterExport(cl=cl, varlist = ls(envir), envir = envir)


  # Estimate model and bootstrap SEs/CIs by method
  if(method=="10CVoneSE"){

    if(perf.metric == "auc" && (class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix")){

      y = as.matrix(data[,1])
      x = as.matrix(data[,-1])

      set.seed(seed)
      fit.roc = glmnet::cv.glmnet(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)

      y.test = as.matrix(test[,1])
      x.test = as.matrix(test[,-1])
      p.fit.test <- stats::predict(fit.roc, s=fit.roc$lambda.1se, newx=x.test, type = "response")
      roc = getROC_AUC(p.fit.test, y.test)
      roc <- list(roc$stack_y, roc$stack_x)

    }

    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])

      set.seed(seed)
      fit = glmnet::cv.glmnet(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)
      coef <- coef(fit, s = "lambda.1se")[-1]

      if(perf.metric == "auc"){

        y.train = as.matrix(data[,1])
        x.train = as.matrix(data[,-1])
        p.fit.train <- stats::predict(fit, s=fit$lambda.1se, newx=x.train, type = "response")
        roc.auc = getROC_AUC(p.fit.train, y.train)
        metric <- roc.auc$auc

      }else{

        metric <- fit$cvm[fit$lambda == fit$lambda.1se]

      }


      # Check if test data is supplied
      if(class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix"){

        # Check compatibility train and test data
        if(dep.var == colnames(test[1]) && ncol(data) == ncol(test)){


          # Test data
          y.test = as.matrix(test[,1])
          x.test = as.matrix(test[,-1])
          p.fit.test <- stats::predict(fit, s=fit$lambda.1se, newx=x.test)
          SSE <- sum((p.fit.test - y.test)^2)

          if(perf.metric == "mse"){

            metric.test <- SSE/nrow(test)

          }else if(perf.metric == "auc"){

            y.test = as.matrix(test[,1])
            x.test = as.matrix(test[,-1])
            p.fit.test <- stats::predict(fit, s=fit$lambda.1se, newx=x.test, type = "response")
            roc.auc = getROC_AUC(p.fit.test, y.test)
            metric.test <- roc.auc$auc

          }else if(perf.metric == "class"){

            y.test = as.matrix(test[,1])
            x.test = as.matrix(test[,-1])
            p.fit.test <- as.numeric(stats::predict(fit, s=fit$lambda.1se, newx=x.test, type = "class"))
            metric.test <- mean(y.test != p.fit.test, na.rm = TRUE)

          }else{

            metric.test <- mean(abs(y.test - p.fit.test), na.rm = TRUE)

          }

          return(c(coef, metric, metric.test))

        }else{

          warning("Training and test data need to have the same outcome and feature variables.")
          stop()

        }

      }else{

        return(c(coef, metric))

      }

    }

    if(seed==0){

      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }else{

      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }


  }else if(method=="10CVmin"){

    if(perf.metric == "auc" && (class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix")){

      y = as.matrix(data[,1])
      x = as.matrix(data[,-1])

      set.seed(seed)
      fit.roc = glmnet::cv.glmnet(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)

      y.test = as.matrix(test[,1])
      x.test = as.matrix(test[,-1])
      p.fit.test <- stats::predict(fit.roc, s=fit.roc$lambda.1se, newx=x.test, type = "response")
      roc = getROC_AUC(p.fit.test, y.test)
      roc <- list(roc$stack_y, roc$stack_x)

    }

    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      set.seed(seed)
      fit = glmnet::cv.glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac, type.measure=perf.metric)
      coef <- coef(fit, s = "lambda.min")[-1]

      if(perf.metric == "auc"){

        y.train = as.matrix(data[,1])
        x.train = as.matrix(data[,-1])
        p.fit.train <- stats::predict(fit, s=fit$lambda.min, newx=x.train, type = "response")
        roc.auc = getROC_AUC(p.fit.train, y.train)
        metric <- roc.auc$auc

      }else{

        metric <- fit$cvm[fit$lambda == fit$lambda.min]

      }


      # Check if test data is supplied
      if(class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix"){

        # Check compatibility train and test data
        if(dep.var == colnames(test[1]) && ncol(data) == ncol(test)){


          # Test data
          y.test = as.matrix(test[,1])
          x.test = as.matrix(test[,-1])
          p.fit.test <- stats::predict(fit, s=fit$lambda.min, newx=x.test)
          SSE <- sum((p.fit.test - y.test)^2)

          if(perf.metric == "mse"){

            metric.test <- SSE/nrow(test)

          }else if(perf.metric == "auc"){

            y.test = as.matrix(test[,1])
            x.test = as.matrix(test[,-1])
            p.fit.test <- stats::predict(fit, s=fit$lambda.1se, newx=x.test, type = "response")
            roc.auc = getROC_AUC(p.fit.test, y.test)
            metric.test <- roc.auc$auc

          }else if(perf.metric == "class"){

            y.test = as.matrix(test[,1])
            x.test = as.matrix(test[,-1])
            p.fit.test <- as.numeric(stats::predict(fit, s=fit$lambda.1se, newx=x.test, type = "class"))
            metric.test <- mean(y.test != p.fit.test, na.rm = TRUE)

          }else{

            metric.test <- mean(abs(y.test - p.fit.test), na.rm = TRUE)

          }

          return(c(coef, metric, metric.test))

        }else{

          warning("Training and test data need to have the same outcome and feature variables.")

        }

      }else{

        return(c(coef, metric))

      }

    }

    if(seed==0){

      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }else{

      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }


  }else{

    boot.glmnet <- function(data, indices, R) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = glmnet::glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac)
      return(coef(fit, s = 0.1)[-1])


    }

    if(seed==0){

      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }else{

      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet, R=r, parallel = "snow", ncpus =cpus, cl=cl)

    }

  }

  # Shut down clusters
  parallel::stopCluster(cl=cl)


  name <- NULL
  coef <- NULL
  SE <- NULL
  KI_low <- NULL
  KI_up <- NULL
  p.val <- NULL
  star <- NULL

  if(type=="norm"){

    type.name="normal"
    bt.s <- 2
    bt.e <- 3

  }else if(type=="perc"){

    type.name="percent"
    bt.s <- 4
    bt.e <- 5

  }else{

    type.name=type
    bt.s <- 4
    bt.e <- 5

  }


  if(method == "none"){

    # Generate output object for method "none"
    iter.loop <- 0

    for(C in c){

      iter.loop <- iter.loop + 1
      coef.nm <- coef.name[iter.loop]
      KI <- boot::boot.ci(results, type=type, conf=conf, index=C)[[type.name]][bt.s:bt.e]

      if(coef.nm %ni% cf.no.shrnkg){

        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- NA
        SE <- c(SE, SE_new)
        p.val_new <- NA
        p.val <- c(p.val, p.val_new)
        KI_low_new <- NA
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- NA
        KI_up <- c(KI_up, KI_up_new)
        star_new <- " "
        star <- c(star, star_new)

      }else{

        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,stats::sd)[C]
        SE <- c(SE, SE_new)
        p.val_new <- mean(abs(results$t0[C]) <= abs(results$t[,C]-mean(results$t[,C])))
        p.val <- c(p.val, p.val_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        star_new <- if(p.val[C] < 0.05 && p.val[C] > 0.01){
          "*"
        }else if(p.val[C] < 0.01 && p.val[C] > 0.001){
          "**"
        }else if(p.val[C] < 0.001){
          "***"
        }else{
          " "
        }
        star <- c(star, star_new)


      }
      output <- list(model_type = model_type,
                     coefficients = name,
                     outcome_var = dep.var,
                     variables_no_shrinkage = cf.no.shrnkg,
                     estimat = coef,
                     standard_error = SE,
                     p.value = p.val,
                     CI_low = KI_low,
                     CI_up = KI_up,
                     star = star,
                     metric_name = "not applied")
    }

  }else{

    # Generate output object for method "10CVmin" and "10CVoneSE"
    if(class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix"){

      metric.KI <- boot::boot.ci(results, type=type, conf=conf, index=max(c)+1)[[type.name]][bt.s:bt.e]
      metric <- results[["t0"]][max(c)+1]

      metric.KI.test <- boot::boot.ci(results, type=type, conf=conf, index=max(c)+2)[[type.name]][bt.s:bt.e]
      metric.test <- results[["t0"]][max(c)+2]
      test.set <- "test data supplied"


    }else{

      metric.KI <- boot::boot.ci(results, type=type, conf=conf, index=max(c)+1)[[type.name]][bt.s:bt.e]
      metric <- results[["t0"]][max(c)+1]
      test.set <- "no test data supplied"

    }

    iter.loop <- 0

    for(C in c){

      iter.loop <- iter.loop + 1
      coef.nm <- coef.name[iter.loop]
      KI <- boot::boot.ci(results, type=type, conf=conf, index=C)[[type.name]][bt.s:bt.e]

      if(coef.nm %ni% cf.no.shrnkg){

        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- NA
        SE <- c(SE, SE_new)
        p.val_new <- NA
        p.val <- c(p.val, p.val_new)
        KI_low_new <- NA
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- NA
        KI_up <- c(KI_up, KI_up_new)
        star_new <- " "
        star <- c(star, star_new)

      }else{

        name_new <- coef.name[iter.loop]
        name <- c(name, name_new)
        coef_new <- results[["t0"]][C]
        coef <- c(coef, coef_new)
        SE_new <- apply(results$t,2,stats::sd)[C]
        SE <- c(SE, SE_new)
        p.val_new <- mean(abs(results$t0[C]) <= abs(results$t[,C]-mean(results$t[,C])))
        p.val <- c(p.val, p.val_new)
        KI_low_new <- KI[1]
        KI_low <- c(KI_low, KI_low_new)
        KI_up_new <- KI[2]
        KI_up <- c(KI_up, KI_up_new)
        star_new <- if(p.val[C] < 0.05 && p.val[C] > 0.01){
          "*"
        }else if(p.val[C] < 0.01 && p.val[C] > 0.001){
          "**"
        }else if(p.val[C] < 0.001){
          "***"
        }else{
          " "
        }
        star <- c(star, star_new)


      }

      if(perf.metric == "auc" && (class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix")){

        output <- list(model_type = model_type,
                       coefficients = name,
                       outcome_var = dep.var,
                       variables_no_shrinkage = cf.no.shrnkg,
                       estimat = coef,
                       standard_error = SE,
                       p.value = p.val,
                       CI_low = KI_low,
                       CI_up = KI_up,
                       star = star,
                       metric_name = perf.metric,
                       metric = metric,
                       metric_CI = metric.KI,
                       test.set = test.set,
                       metric.test = metric.test,
                       metric_CI.test = metric.KI.test,
                       roc_y = roc[1],
                       roc_x = roc[2])

      }else if(class(test)== "data.frame" || class(test)== "tbl_df" || class(test)== "matrix"){

        output <- list(model_type = model_type,
                       coefficients = name,
                       outcome_var = dep.var,
                       variables_no_shrinkage = cf.no.shrnkg,
                       estimat = coef,
                       standard_error = SE,
                       p.value = p.val,
                       CI_low = KI_low,
                       CI_up = KI_up,
                       star = star,
                       metric_name = perf.metric,
                       metric = metric,
                       metric_CI = metric.KI,
                       test.set = test.set,
                       metric.test = metric.test,
                       metric_CI.test = metric.KI.test)

      }else{

        output <- list(model_type = model_type,
                       coefficients = name,
                       outcome_var = dep.var,
                       variables_no_shrinkage = cf.no.shrnkg,
                       test.set = test.set,
                       estimat = coef,
                       standard_error = SE,
                       p.value = p.val,
                       CI_low = KI_low,
                       CI_up = KI_up,
                       star = star,
                       metric_name = perf.metric,
                       metric = metric,
                       metric_CI = metric.KI)

      }

    }

  }

  class(output) <- c("glmnetSE")

  invisible(output)

}

#####################
# summary function #
#####################

#' @title Summary Function for a fitted glmnetSE Objects
#' @description Print the coefficients with standard errors, confidence intervals, and p-values of a \code{\link{glmnetSE}} model. The inference statistics are only available for the coefficients without shrinkage applied. They would be biased otherwise. Only if cross-fold validation is used in the \code{glmnetSE} model, the selected performance metric is displayed. If test data is supplied the performance metric on the train as test data is displayed.
#' @param object A model of the class \code{glmnetSE}.
#' @param ... Additional arguments affecting the summary produced.
#' @return The output of a \code{glmnetSE} object and the performance metric if cross-fold validation is used.
#' @keywords glmnetSE summary results output
#' @examples
#' \dontrun{
#' # Estimate model
#'
#' glmnetSE.model <- glmnetSE(data=swiss,cf.no.shrnkg = c("Education"))
#'
#'
#' # Display model output with summary
#'
#' summary(glmnetSE.model)
#' }
#' @export

summary.glmnetSE <- function(object, ...){

  # Warning and stop if object not of the class 'glmnetSE' is used
  if(class(object)!= "glmnetSE"){
    warning("A object of the class 'glmnetSE' should be used")
    stop()
  }

  if(object$metric_name == "not applied"){

    coef <- data.frame(cbind("Coefficients" = object$coefficients,
                             "Estimates" = round(object$estimat,4),
                             "Std.Error" = round(object$standard_error,4),
                             "CI.low" = round(object$CI_low,4),
                             "CI.up" = round(object$CI_up,4),
                             "p-value" = round(object$p.value,4),
                             "Sig." = object$star))



  }else{

    coef <- data.frame(cbind("Coefficients" = object$coefficients,
                             "Estimates" = round(object$estimat,4),
                             "Std.Error" = round(object$standard_error,4),
                             "CI.low" = round(object$CI_low,4),
                             "CI.up" = round(object$CI_up,4),
                             "p-value" = round(object$p.value,4),
                             "Sig." = object$star))

    if(object$test.set == "test data supplied"){

      metric <- data.frame(cbind("Metric" = c(object$metric_name),
                                 "Value" = c(round(object$metric,2)),
                                 "CI.low" = c(round(object$metric_CI[1],2)),
                                 "CI.up" = c(round(ifelse(object$metric_name == "auc" && object$metric_CI[2] > 1,1, object$metric_CI[2]),2))))

      metric.test <- data.frame(cbind("Metric" = c(object$metric_name),
                                      "Value" = c(round(object$metric.test,2)),
                                      "CI.low" = c(round(object$metric_CI.test[1],2)),
                                      "CI.up" = c(round(ifelse(object$metric_name == "auc" && object$metric_CI[2] > 1,1, object$metric_CI[2]),2))))


    }else{

      metric <- data.frame(cbind("Metric" = c(object$metric_name),
                                 "Value" = c(round(object$metric,2)),
                                 "CI.low" = c(round(object$metric_CI[1],2)),
                                 "CI.up" =  c(round(ifelse(object$metric_name == "auc" && object$metric_CI[2] > 1,1, object$metric_CI[2]),2))))

    }

  }

  cat(paste0(object$model_type), "model:\n")
  cat("Outcome variable:", paste0(object$outcome_var), "\n")
  cat("Variables without shrinkage:", paste0(object$variables_no_shrinkage), "\n")
  cat("\n\n")
  print(coef, row.names = FALSE)
  cat("\n\n")

  if(object$metric_name == "not applied"){

    cat("Performance Metric: \n")
    cat(c("Metric only applicable with method '10CVmin' or '10CVoneSE'"), fill = getOption("width"))

  }else{

    if(object$test.set == "test data supplied"){

      cat("Performance metric train data: \n")
      print(metric, row.names = FALSE)
      cat("\n\n")
      cat("Performance metric test data: \n")
      print(metric.test, row.names = FALSE)

    }else{

      cat("Performance Metric: \n")
      print(metric, row.names = FALSE)

    }

  }

}

#####################
# ROC plot function #
#####################

#' @title Plot ROC Curve of a fitted glmnetSE Model on Test Data
#' @description Plot the ROC curve of a fitted  model \code{\link{glmnetSE}} (family \code{binomial} and performance metric \code{auc}) on supplied test data.
#' @param x A model of the class \code{glmnetSE} of family \code{binomial} and performance metric \code{auc} for which the ROC curve should be plotted.
#' @param ... Additional arguments affecting the plot produced.
#' @return The ROC curve of a \code{glmnetSE} object.
#' @keywords glmnetSE plot ROC
#' @examples
#' \dontrun{
#' #' # Generate dichotom variable
#'
#' swiss$Fertility <- ifelse(swiss$Fertility >= median(swiss$Fertility), 1, 0)
#'
#'
#' # Estimate model
#'
#' glmnetSE.model <- glmnetSE(data=swiss.train, cf.no.shrnkg = c("Education"),
#' method = "10CVoneSE", test = swiss.test, family = "binomial", perf.metric = "auc")
#'
#'
#' # Plot ROC curve of the fitted model on swiss.test data
#'
#' plot(glmnetSE.model)
#' }
#' @export


plot.glmnetSE <- function(x, ...){

  # Warning and stop if object not of the class 'glmnetSE' is used
  if(class(x)!= "glmnetSE"){
    warning("A object of the class 'glmnetSE' should be used")
    stop()
  }

  # Warning and stop if no test data is supplied.
  if(x$test.set!= "test data supplied"){
    warning("The ROC curve can only be plotted if test data is supplied to 'glmnetSE'")
    stop()
  }

  graphics::plot(x[["roc_x"]][[1]], x[["roc_y"]][[1]], type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC")
  graphics::axis(1, seq(0.0,1.0,0.1))
  graphics::axis(2, seq(0.0,1.0,0.1))
  graphics::abline(0, 1, col="black", lty=3)

}

