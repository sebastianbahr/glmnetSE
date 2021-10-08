#' @title glmnetSE: Add bootstrap SE to glmnet for selected coefficients (no shrinkage)
#' @description glmnet or cv.glmnet with bootstrap standard errors for selected coefficients with no shrinkage applied for them.
#' @param data A data frame, tibble or matrix object with the outcome variable in the first column and the feature variables in the following columns.
#' @param cf.no.shrnkg A character string of the coefficients that will be interpreted, the inference statistic is of interest and therefore no shrinkage will be applied.
#' @param alpha Alpha value 1 for LASSO and 0 for Ridge.
#' @param method A character string defining if no 10-fold cross validation "none", 10-fold cross validation and selecting the lambda at which the MSE is within one standard error of the smallest MSE  "10CVoneSE", or 10-fold cross validation and selecting the lambda at which the smallest MSE is achieved "10CVmin" - default is "none.
#' @param r Number of bootstrap repetitions -default is 100.
#' @param nlambda Number of tested lambda values - default is 100.
#' @param seed Seed set for the bootstrap sampling - default 0 which means no seed set.
#' @param family A character string representing the used model family either "gaussian" or "binomial" - default is "gaussian".
#' @param type A character string indicating the type of calculated bootstrap intervals. It can be "norm", "basic", "stud", "perc", "bca". For more information check the "boot.ci" package - default is "basic".
#' @param conf Indicates the confidence interval level - default is 0.95.
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

glmnetSE <- function(data, cf.no.shrnkg, alpha=1, method="none", r=100, nlambda=100, seed=0, family="gaussian", type="basic", conf=0.95){


  coef.name <- colnames(data)[-1]
  c <- which(names(data) %in% coef.name)-1
  c <- 1:(ncol(data)-1)
  no.shrink <- which(names(data) %in% cf.no.shrnkg)-1

  p.fac <- rep(1, ncol(data)-1)
  p.fac[no.shrink] <- 0


  "%ni%" <- Negate("%in%")

  if(method=="10CVoneSE"){
    boot.glmnet <- function(data, indices, type, conf) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = glmnet::cv.glmnet(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac)
      return(coef(fit, s = "lambda.1se")[-1])
    }

    if(seed==0){
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }else{
      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }


  }else if(method=="10CVmin"){
    boot.glmnet <- function(data, indices, type, conf) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = glmnet::cv.glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac)
      return(coef(fit, s = "lambda.min")[-1])
    }

    if(seed==0){
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }else{
      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }


  }else{
    boot.glmnet <- function(data, indices, type, conf) {
      d = data[indices,]
      y = as.matrix(d[,1])
      x = as.matrix(d[,-1])
      fit = glmnet::glmnet(x=x, y=y,alpha=alpha, family=family, nlambda=nlambda, penalty.factor=p.fac)
      return(coef(fit, s = 0.1)[-1])
    }

    if(seed==0){
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }else{
      set.seed(seed)
      results = boot::boot(data=data, statistic=boot.glmnet,
                           R=r, type=type, conf=conf)
    }

  }


  output <- data.frame(name = character(), coef = numeric(), KI_low = numeric(), KI_up = numeric(), sig = character())
  iter.loop <- 0
  for(C in c){
    iter.loop <- iter.loop + 1
    coef.nm <- coef.name[iter.loop]
    KI <- boot::boot.ci(results, type="basic", index=C)[["basic"]][4:5]
    if(coef.nm %ni% cf.no.shrnkg){
      name <- coef.name[iter.loop]
      coef <- results[["t0"]][C]
      KI_low <- NA
      KI_up <- NA
      sig <- c("")
    }else if((KI[1] > 0 & KI[2] > 0) | (KI[1] < 0 & KI[2] < 0)){
      name <- coef.name[iter.loop]
      coef <- results[["t0"]][C]
      KI_low <- KI[1]
      KI_up <- KI[2]
      sig <- c("*")
    }else{
      name <- coef.name[iter.loop]
      coef <- results[["t0"]][C]
      KI_low <- KI[1]
      KI_up <- KI[2]
      sig <- c("")
    }
    "%>%" = magrittr::"%>%"
    output <- output %>% dplyr::add_row(name, coef, KI_low, KI_up, sig)
  }
  return(output)


}

