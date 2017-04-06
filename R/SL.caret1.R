#' SL.caret1 
#' 
#' A modification of the SL.caret function from the SuperLearner package that
#' suppresses some output when method = "gbm" or "nnet". See \code{?SL.caret}
#' for more information.
#' 
#' @param Y Training outcomes
#' @param X Training predictors 
#' @param newX Test set predictors 
#' @param obsWeights Weights for the observations
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param family Character indicating family argument (ignored)
#' @param ... Other arguments (not currently used)
#' 
#' @importFrom caret train
#' @importFrom stats predict
#' 
#' @export

SL.caret1 <- function (Y, X, newX, family, obsWeights, method = "rf", tuneLength = 3, 
                       trControl = caret::trainControl(method = "cv", number = ifelse(length(Y)<=100, 20, 10), verboseIter = FALSE), 
                       ...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,verbose=FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    if(method=="gbm"){
      suppressWarnings(
        # pass verbose==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl, verbose = FALSE)
      )
    }else if(method=="nnet"){
      suppressWarnings(
        # pass trace==FALSE directly to train (verboseIter doesn't 
        # suppress output)
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,trace=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}

