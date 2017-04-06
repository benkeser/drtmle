#' SL.rpart.caret1 
#' 
#' Uses SL.caret1 to train regression tree. 
#' 
#' @param method Character describing what algorithm to train
#' @param tuneLength The number of tuning parameter combinations
#' @param trControl Passed to \code{caret::train}
#' @param ... Other arguments passed to SL.caret1
#' 
#' @export

SL.rpart.caret1 <- function(...,method="rpart",tuneLength = 8, trControl = caret::trainControl(method = "cv", number = 5)){
  SL.caret1(...,method=method,tuneLength=tuneLength, trControl = trControl)
}