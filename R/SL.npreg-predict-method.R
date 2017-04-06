#' predict.SL.npreg
#' 
#' Method for predicting SL.npreg
#' 
#' @param object An object of class "SL.npreg"
#' @param newdata The new data used to obtain predictions 
#' @param ... Other arguments passed to predict
#' @importsFrom stats predict
#' @export
predict.SL.npreg <- function(object, newdata, ...){
  pred <- stats::predict(object=object$object, newdata=newdata)
  pred
}