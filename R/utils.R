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


#' Print the output of ci.drtmle
#' @export
#' @param x An object of class ci.drtmle
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print ci.drtmle	

print.ci.drtmle <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}

#' Print the output of ci.islptw
#' @export
#' @param x An object of class ci.islptw
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print ci.islptw 

print.ci.islptw <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}

#' Print the output of wald_test.drtmle
#' @export
#' @param x An object of class wald_test.drtmle
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print wald_test.drtmle

print.wald_test.drtmle <- function(x,digits = 3,...){
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
}

#' Print the output of wald_test.islptw
#' @export
#' @param x An object of class wald_test.islptw
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print wald_test.islptw 

print.wald_test.islptw <- function(x,digits = 3,...){
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
}