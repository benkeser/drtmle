#' Super learner wrapper for kernel regression
#' 
#' Kernel regression based on the \href{https://CRAN.R-project.org/package=np}{np}
#' package. Uses leave-one-out cross-validation to fit a kernel regression. 
#' See \code{?npreg} for more details.
#' 
#' @param Y A vector of outcomes.
#' @param X A matrix or data.frame of training data predictors. 
#' @param family Not used by the function directly, but ensures compatibility with \code{SuperLearner}.
#' @param newX A test set of predictors.
#' @param obsWeights Not used by the function directly, but ensures compatibility with \code{SuperLearner}.
#' @param rangeThresh If the the range of the outcomes is smaller than this number, the method
#' returns the empirical average of the outcomes. Used for computational expediency and stability.  
#' @param ... Other arguments (not currently used).
#' 
#' @importFrom np npregbw npreg
#' @importFrom stats predict as.formula
#' 
#' @export 
#' 
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' X <- data.frame(X1 = rnorm(n))
#' Y <- X$X1 + rnorm(n)
#' # fit npreg
#' fit <- SL.npreg(Y = Y, X = X, newX = X)

SL.npreg <- function (Y, X, newX, family = gaussian(), 
                      obsWeights = rep(1, length(Y)), 
                      rangeThresh=1e-7, ...) 
{
  # turn off annoying messages
  options(np.messages = FALSE)
  # check range threshold
  if(abs(diff(range(Y))) <= rangeThresh){
    thisMod <- glm(Y ~ 1, data = X)
  }else{
    # bandwidth selection
    bw <- np::npregbw(stats::as.formula(
            paste("Y ~", paste(names(X), collapse = "+"))), data = X,
          ftol=0.01, tol=0.01, remin=FALSE)
  
    # fit the kernel regression
    thisMod <- np::npreg(bw)
  }
  
  pred <- stats::predict(thisMod, newdata = newX)
  fit <- list(object = thisMod)
  class(fit) <- "SL.npreg"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Predict method for SL.npreg
#' 
#' Method for predicting SL.npreg objects. 
#' 
#' @param object An object of class \code{"SL.npreg"}.
#' @param newdata The new data used to obtain predictions.
#' @param ... Other arguments passed to predict.
#' 
#' @importFrom stats predict
#' @method predict SL.npreg
#' @export
#' @examples
#' # simulate data
#' set.seed(1234)
#' n <- 100
#' X <- data.frame(X1 = rnorm(n))
#' Y <- X$X1 + rnorm(n)
#' # fit npreg
#' fit <- SL.npreg(Y = Y, X = X, newX = X)
#' # predict on fit
#' newX <- data.frame(X1 = c(-1,0,1))
#' pred <- predict(fit$fit, newdata = newX)

predict.SL.npreg <- function(object, newdata, ...){
  pred <- stats::predict(object = object$object, newdata = newdata)
  pred
}
