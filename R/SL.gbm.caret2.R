#' SL.gbm.caret2
#' 
#' Uses SL.caret 1 to train a gbm with 10 choices of tuning parameters and 5-fold CV
#' 
#' @param trControl Passed to \code{caret::train}
#' @param tuneLength The number of tuning parameter combinations
#' @param method Character describing what algorithm to train
#' @param ... Other arguments passed to \code{SL.caret1}
#' 
#' @export

SL.gbm.caret2 <- function (..., method = "gbm", tuneLength = 20, trControl = caret::trainControl(method = "cv", number = 5)) 
{
    SL.caret1(..., method = method, tuneLength = tuneLength, trControl = trControl)
}
