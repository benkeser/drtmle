#--------------------------------------------------
# Print methods
#--------------------------------------------------
#' Print the output of a \code{"drtmle"} object.
#' 
#' @param x A \code{"drtmle"} object
#' @param ... Other arguments (not used)
#' @export
#' @method print drtmle

print.drtmle <- function(x,...){
  tmp <- list(est = cbind(x$drtmle$est),
  	          cov = x$drtmle$cov)
  row.names(tmp$est) <- x$a_0
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if(length(x$a_0) <= 4){
  	print(tmp)
  }else{
  	tmp$cov <- diag(tmp$cov)
  }
}

#' Print the output of a \code{"islptw"} object.
#' 
#' @param x A \code{"islptw"} object.
#' @param ... Other arguments (not used)
#' @export
#' @method print islptw

print.islptw <- function(x,...){
  tmp <- list(est = cbind(x$islptw_tmle$est),
  	          cov = x$islptw_tmle$cov)
  row.names(tmp$est) <- x$a_0
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if(length(x$a_0) <= 4){
  	print(tmp)
  }else{
  	tmp$cov <- diag(tmp$cov)
  }
}


#' Print the output of confint.drtmle
#' @export
#' @param x An object of class drconfint
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print confint.drtmle	

print.confint.drtmle <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}

#' Print the output of confint.islptw
#' @export
#' @param x An object of class drconfint
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print confint.islptw 

print.confint.islptw <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}