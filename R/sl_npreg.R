
#' Super learner wrapper for kernel regression
#' 
#' Kernel regression based on the \href{https://cran.r-project.org/web/packages/np/}{np}
#' package. Uses leave-one-out cross-validation to fit a kernel regression. 
#' See \code{?npreg} for more details.
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param family Not used by the function, but needed for compatibility with \code{SuperLearner}
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param rangeThresh For stability: if the the range of the outcomes is smaller than this number, the method
#' returns the empirical average of the outcomes. 
#' @param ... Other arguments (not currently used)
#' 
#' @importFrom np npregbw npreg
#' @importFrom stats predict as.formula
#' 
#' @export 

SL.npreg <- function (Y, X, newX, family, obsWeights, 
                      rangeThresh=1e-7, ...) 
{
  if(abs(diff(range(Y))) <= rangeThresh){
    thisMod <- glm(Y ~ 1, data=X)
  }else{
  bw <- np::npregbw(stats::as.formula(paste("Y ~", paste(names(X),collapse="+"))), data=X,
                ftol=0.01, tol=0.01, remin=FALSE)
  
  # fit the kernel regression
  thisMod <- np::npreg(bw)
  }
  
  pred <- stats::predict(thisMod, newdata=newX)
  fit <- list(object = thisMod)
  class(fit) <- "SL.npreg"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Predict method for SL.npreg
#' 
#' Method for predicting SL.npreg objects. 
#' 
#' @param object An object of class \code{"SL.npreg"}
#' @param newdata The new data used to obtain predictions 
#' @param ... Other arguments passed to predict
#' 
#' @importFrom stats predict
#' @method predict SL.npreg
#' @export
predict.SL.npreg <- function(object, newdata, ...){
  pred <- stats::predict(object=object$object, newdata=newdata)
  pred
}
