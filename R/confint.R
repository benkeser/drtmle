#' Confidence intervals for drtmle objects
#' 
#' Estimate doubly-robust confidence intervals for
#' an object of class \code{"drtmle"}
#' 
#' @param object An object of class \code{"drtmle"}
#' @param parm A vector indicating for which estimators to return a 
#' confidence interval. Possible estimators include the TMLE with doubly robust
#' inference (\code{"drtmle"}, recommended), the AIPTW with additional correction 
#' for misspecification (\code{"aiptw_c"}, not recommended), the standard TMLE 
#' (\code{"tmle"}, recommended only for comparison to "drtmle"), the standard 
#' AIPTW (\code{"aiptw"}, recommended only for comparison to "drtmle").  
#' @param level The level of the confidence interval
#' @param contrast This option specifies what parameter to return confidence intervals for.
#' If \code{contrast=NULL}, then compute confidence intervals for the covariate-adjusted
#' marginal means. \code{contrast} can also be input as a numeric vector of ones, negative ones,
#' and zeros to define linear combinations of the various means (e.g., to estimate an average
#' treatment effect, see examples). \code{contrast} can also be a list with named functions
#' \code{f}, \code{f_inv}, \code{h}, and \code{h_grad}. The first two functions should take
#' as input argument \code{eff}. Respectively, these specify which transformation 
#' of the effect measure to compute the confidence interval for and the inverse
#' transformation to put the confidence interval back on the original scale. The function \code{h}
#' defines the contrast to be estimated and should take as input \code{parm}, a vector
#' of the same length as \code{object$a_0}, and output the desired contrast. The function
#' \code{h_grad} is the gradient of the function \code{h}. See examples. 
#' @param ... Other options (not currently used)
#' @importFrom stats qnorm
#' @export
#' @method confint drtmle
#' @return An object of class \code{"confint.drtmle"} with point estimates and
#' confidence intervals of the specified level. 
#' 
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # simulate data
#' set.seed(123456)
#' n <- 200
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#' # fit drtmle
#' fit1 <- drtmle(W = W, A = A, Y = Y, a_0 = c(1,0),
#'             family=binomial(),
#'             stratify=FALSE,
#'             SL_Q=c("SL.glm","SL.mean","SL.step"),
#'             SL_g=c("SL.glm","SL.mean","SL.step"),
#'             SL_Qr="SL.npreg",
#'             SL_gr="SL.npreg")
#' 
#' # get confidence intervals for each mean
#' ci_mean <- confint(fit1)
#' 
#' # get confidence intervals for ATE
#' ci_ATE <- confint(fit1, contrast = c(1,-1))
#' 
#' # get confidence intervals for risk ratio by
#' # computing CI on log scale and back-transforming
#' myContrast <- list(f = function(eff){ log(eff) },
#'                    f_inv = function(eff){ exp(eff) },
#'                    h = function(parm){ parm[1]/parm[2] },
#'                    h_grad =  function(parm){ c(1/parm[1],-1/parm[2]) })
#' ci_RR <- confint(fit1, contrast = myContrast)


confint.drtmle <- function(object, parm = c("drtmle"), level = 0.95,
                      contrast = NULL,...){
	if(class(object) != "drtmle"){
		stop("drconfint only works with drtmle objects")
	}
	out <- vector(mode = "list", length = length(parm))
	names(out) <- parm
	# if no contrast then return an CI for each 
	# covariate-adjusted mean
	if(is.null(contrast)){
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 3)
			for(j in seq_along(object$a_0)){
				out[[i]][j,] <-
				  rep(object[[parm[i]]]$est[j],3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				  	rep(sqrt(object[[parm[i]]]$cov[j,j]),3)
			}
			row.names(out[[i]]) <- object$a_0
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}else if(is.numeric(contrast)){
		# check that contrast came in correctly
		if(length(contrast) != length(object$a_0)){
			stop("length of contrast vector not equal to length of a_0")
		}
		if(!all(contrast %in% c(-1,1,0))){
			stop("contrast should only be -1, 1, or 0")
		}
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
			g <- matrix(contrast, nrow = length(contrast))
			p <- matrix(object[[parm[i]]]$est, nrow = length(contrast))
			v <- object[[parm[i]]]$cov
			thisC <- t(g)%*%p
			thisSe <- sqrt(t(g)%*%v%*%g)
			out[[i]][1,] <- rep(thisC,3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				rep(thisSe,3)
			# E[Y(a_0[1])] - E[Y(a_0[2])]
			indMinus <- which(contrast == -1)
			indPlus <- which(contrast == 1)
			plusTerms <- minusTerms <- ""
			if(length(indMinus) > 0){
				minusTerms <- paste0("E[Y(",object$a_0[indMinus],")]")
			}
			if(length(indPlus) > 0){
				plusTerms <- paste0("E[Y(",object$a_0[indPlus],")]")
			}
			thisName <- paste0(paste0(plusTerms,collapse = "+"),
			                   ifelse(length(indMinus)>0,"-",""),
			                   paste0(minusTerms,collapse="-"))

			row.names(out[[i]]) <- thisName
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}else if(is.list(contrast)){
		if(!all(c("f","f_inv","h","h_grad") %in% names(contrast))){
			stop("some function missing in contrast. see ?drconfint for help.")
		}
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
			thisC <- do.call(contrast$h,args=list(parm = object[[parm[i]]]$est))
			f_thisC <- do.call(contrast$f,args=list(eff = thisC))

			grad <- matrix(do.call(contrast$h_grad,
			                       args=list(parm = object[[parm[i]]]$est)
			               ),nrow = length(object$a_0))
			v <- object[[parm[i]]]$cov
			thisSe <- sqrt(t(grad)%*%v%*%grad)
			transformCI <- rep(f_thisC,3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				rep(thisSe,3)
			out[[i]][1,] <- do.call(contrast$f_inv, args = list(eff = transformCI))
			row.names(out[[i]]) <- c("user contrast")
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}
	class(out) <- "confint.drtmle"
	return(out)
}

#' Confidence intervals for islptw objects
#' 
#' Estimate confidence intervals for objects of class \code{"islptw"}
#' 
#' @param object An object of class \code{"islptw"}
#' @param parm A vector indicating for which estimators to return a 
#' confidence interval. Possible estimators include the TMLE IPTW 
#' (\code{"islptw_tmle"}, recommended), the one-step IPTW 
#' (\code{"islptw_os"}, not recommended), the standard IPTW 
#' (\code{"iptw"}, recommended only for comparison to the other two estimators).
#' @param level The level of the confidence interval
#' @param contrast This option specifies what parameter to return confidence intervals for.
#' If \code{contrast=NULL}, then compute confidence intervals for the covariate-adjusted
#' marginal means. \code{contrast} can also be input as a numeric vector of ones, negative ones,
#' and zeros to define linear combinations of the various means (e.g., to estimate an average
#' treatment effect, see examples). \code{contrast} can also be a list with named functions
#' \code{f}, \code{f_inv}, \code{h}, and \code{h_grad}. The first two functions should take
#' as input argument \code{eff}. Respectively, these specify which transformation 
#' of the effect measure to compute the confidence interval for and the inverse
#' transformation to put the confidence interval back on the original scale. The function \code{h}
#' defines the contrast to be estimated and should take as input \code{parm}, a vector
#' of the same length as \code{object$est}, and output the desired contrast. The function
#' \code{h_grad} is the gradient of the function \code{h}. See examples. 
#' @param ... Other options (not currently used)

#' @export
#' @method confint islptw
#' @return An object of class \code{"confint.islptw"} with point estimates and
#' confidence intervals of the specified level. 
#' 
#' @examples
#' # load super learner
#' library(SuperLearner)
#' # fit islptw
#' set.seed(123456)
#' n <- 200
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#'
#' fit1 <- islptw(W = W, A = A, Y = Y, a_0 = c(1,0),
#'                SL_g=c("SL.glm","SL.mean","SL.step"),
#'                SL_Qr="SL.glm")
#' 
#' # get confidence intervals for each mean
#' ci_mean <- confint(fit1)
#' 
#' # get confidence intervals for ATE
#' ci_ATE <- confint(fit1, contrast = c(1,-1))
#' 
#' # get confidence intervals for risk ratio 
#' # by inputting own contrast function
#' # this computes CI on log scale and back transforms
#' myContrast <- list(f = function(eff){ log(eff) },
#'                    f_inv = function(eff){ exp(eff) },
#'                    h = function(parm){ parm[1]/parm[2] },
#'                    h_grad =  function(parm){ c(1/parm[1],-1/parm[2]) })
#' ci_RR <- confint(fit1, contrast = myContrast)


confint.islptw <- function(object, parm = c("islptw_tmle"), level = 0.95,
                      contrast = NULL,...){
	if(any(parm=="iptw")){
		stop("Theory does not support inference for naive IPTW with super learner.")
	}
	if(class(object) != "islptw"){
		stop("drconfint only works with islptw objects")
	}
	out <- vector(mode = "list", length = length(parm))
	names(out) <- parm
	# if no contrast then return an CI for each 
	# covariate-adjusted mean
	if(is.null(contrast)){
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 3)
			for(j in seq_along(object$a_0)){
				out[[i]][j,] <-
				  rep(object[[parm[i]]]$est[j],3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				  	rep(sqrt(object[[parm[i]]]$cov[j,j]),3)
			}
			row.names(out[[i]]) <- object$a_0
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}else if(is.numeric(contrast)){
		# check that contrast came in correctly
		if(length(contrast) != length(object$a_0)){
			stop("length of contrast vector not equal to length of a_0")
		}
		if(!all(contrast %in% c(-1,1,0))){
			stop("contrast should only be -1, 1, or 0")
		}
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
			g <- matrix(contrast, nrow = length(contrast))
			p <- matrix(object[[parm[i]]]$est, nrow = length(contrast))
			v <- object[[parm[i]]]$cov
			thisC <- t(g)%*%p
			thisSe <- sqrt(t(g)%*%v%*%g)
			out[[i]][1,] <- rep(thisC,3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				rep(thisSe,3)
			# E[Y(a_0[1])] - E[Y(a_0[2])]
			indMinus <- which(contrast == -1)
			indPlus <- which(contrast == 1)
			plusTerms <- minusTerms <- ""
			if(length(indMinus) > 0){
				minusTerms <- paste0("E[Y(",object$a_0[indMinus],")]")
			}
			if(length(indPlus) > 0){
				plusTerms <- paste0("E[Y(",object$a_0[indPlus],")]")
			}
			thisName <- paste0(paste0(plusTerms,collapse = "+"),
			                   ifelse(length(indMinus)>0,"-",""),
			                   paste0(minusTerms,collapse="-"))

			row.names(out[[i]]) <- thisName
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}else if(is.list(contrast)){
		if(!all(c("f","f_inv","h","h_grad") %in% names(contrast))){
			stop("some function missing in contrast.")
		}
		for(i in seq_along(parm)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 3)
			thisC <- do.call(contrast$h,args=list(parm = object[[parm[i]]]$est))
			f_thisC <- do.call(contrast$f,args=list(eff = thisC))

			grad <- matrix(do.call(contrast$h_grad,
			                       args=list(parm = object[[parm[i]]]$est)
			               ),nrow = length(object$a_0))
			v <- object[[parm[i]]]$cov
			thisSe <- sqrt(t(grad)%*%v%*%grad)
			transformCI <- rep(f_thisC,3) + stats::qnorm(c(0.5,(1-level)/2,(1+level)/2))*
				rep(thisSe,3)
			out[[i]][1,] <- do.call(contrast$f_inv, args = list(eff = transformCI))
			row.names(out[[i]]) <- c("user contrast")
			colnames(out[[i]]) <- c("est","cil","ciu")
		}
	}
	class(out) <- "confint.islptw"
	return(out)
}