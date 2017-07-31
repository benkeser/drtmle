#' Wald tests for drtmle and islptw objects
#' @param ... Arguments to be passed to method
#' @export
wald_test <- function(...){
	UseMethod("wald_test")
}

#' Wald tests for drtmle objects
#' 
#' @param object An object of class \code{"drtmle"}
#' @param est A vector indicating for which estimators to return a 
#' confidence interval. Possible estimators include the TMLE with doubly robust
#' inference (\code{"drtmle"}, recommended), the AIPTW with additional correction 
#' for misspecification (\code{"aiptw_c"}, not recommended), the standard TMLE 
#' (\code{"tmle"}, recommended only for comparison to "drtmle"), the standard 
#' AIPTW (\code{"aiptw"}, recommended only for comparison to "drtmle"), and 
#' G-computation (\code{"gcomp"}, not recommended).
#' @param null The null hypothesis value. 
#' @param contrast This option specifies what parameter to return confidence intervals for.
#' If \code{contrast=NULL}, then test the null hypothesis that the covariate-adjusted
#' marginal means equal the value(s) specified in \code{null}. 
#' \code{contrast} can also be a numeric vector of ones, negative ones,
#' and zeros to define linear combinations of the various means (e.g., to estimate an average
#' treatment effect, see examples). In this case, we test the null hypothesis that the 
#' linear combination of means equals the value specified in \code{null}. 
#' \code{contrast} can also be a list with named functions
#' \code{f}, \code{h}, and \code{h_grad}. The function \code{f} takes
#' as input argument \code{eff} and specifies which transformation 
#' of the effect measure to test. The function \code{h} defines the contrast 
#' to be estimated and should take as input \code{est}, a vector
#' of the same length as \code{object$a_0}, and output the desired contrast. The function
#' \code{h_grad} is the gradient of the function \code{h(f())}. The function
#' computes a test of the null hypothesis that \code{h(f(object$est)) = null}. 
#' See examples. 
#' @param ... Other options (not currently used)
#' 
#' @importFrom stats pnorm
#' @export
#' @method wald_test drtmle
#' @return An object of class \code{"ci.drtmle"} with point estimates and
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
#' # get hypothesis test that each mean = 0.5
#' test_mean <- wald_test(fit1, null = 0.5)
#' 
#' # get test that ATE = 0
#' test_ATE <- wald_test(fit1, null = 0, contrast = c(1,-1))
#' 
#' # get test that risk ratio = 1, computing test on log scale
#' myContrast <- list(f = function(eff){ log(eff) },
#'                    f_inv = function(eff){ exp(eff) },
#'                    h = function(est){ est[1]/est[2] },
#'                    h_grad =  function(est){ c(1/est[1],-1/est[2]) })
#' ci_RR <- ci(fit1, contrast = myContrast, null = 1)

wald_test.drtmle <- function(object, est = c("drtmle"), null = 0,
                             contrast = NULL, ...){
	if(class(object) != "drtmle"){
		stop("wald_test only works with drtmle objects")
	}
	out <- vector(mode = "list", length = length(est))
	names(out) <- est
	# if no contrast then return an CI for each 
	# covariate-adjusted mean
	if(is.null(contrast)){
		# make sure null is input properly
		if(length(null) == 1){
			null <- rep(null, length(object$a_0))
		}else if(length(null) > length(object$a_0)){
			stop("length of null should be one or same length as object$a_0")
		}

		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 2)
			for(j in seq_along(object$a_0)){
				zstat <- (object[[est[i]]]$est[j]-null[j])/sqrt(object[[est[i]]]$cov[j,j])
				pval <- 2*stats::pnorm(-abs(zstat))
				out[[i]][j,] <- c(zstat,pval)
			}
			row.names(out[[i]]) <- paste0(
              "H0: E[Y(",object$a_0,")]=",null
          	)
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}else if(is.numeric(contrast)){
		# check that contrast came in correctly
		if(length(contrast) != length(object$a_0)){
			stop("length of contrast vector not equal to length of a_0")
		}
		if(!all(contrast %in% c(-1,1,0))){
			stop("contrast should only be -1, 1, or 0")
		}
		if(length(null) > 1){
			stop("null length should be 1 if contrast is a linear combination")
		}
		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 2)
			g <- matrix(contrast, nrow = length(contrast))
			p <- matrix(object[[est[i]]]$est, nrow = length(contrast))
			v <- object[[est[i]]]$cov
			thisC <- t(g)%*%p
			thisSe <- sqrt(t(g)%*%v%*%g)
			zstat <- (thisC - null) / thisSe
			pval <- 2*stats::pnorm(-abs(zstat))
			out[[i]][1,] <- c(zstat, pval)
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
			thisName <- paste0("H0:",paste0(plusTerms,collapse = "+"),
			                   ifelse(length(indMinus)>0,"-",""),
			                   paste0(minusTerms,collapse="-"),"=",null)

			row.names(out[[i]]) <- thisName
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}else if(is.list(contrast)){
		if(!all(c("f","h","h_grad") %in% names(contrast))){
			stop("some function missing in contrast. see ?wald_test for help.")
		}
		if(length(null) > 1){
			stop("null length should be 1 if contrast is user-specified functions")
		}
		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 2)
			thisC <- do.call(contrast$h,args=list(est = object[[est[i]]]$est))
			f_thisC <- do.call(contrast$f,args=list(eff = thisC))

			grad <- matrix(do.call(contrast$h_grad,
			                       args=list(est = object[[est[i]]]$est)
			               ),nrow = length(object$a_0))
			v <- object[[est[i]]]$cov
			thisSe <- sqrt(t(grad)%*%v%*%grad)
			f_null <- do.call(contrast$f, args = list(eff = null))
			zstat <- (f_thisC - f_null)/ thisSe
			pval <- 2*stats::pnorm(-abs(zstat))
			out[[i]][1,] <- c(zstat,pval)
			row.names(out[[i]]) <- paste0("H0: user contrast = ",null)
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}
	class(out) <- "wald_test.drtmle"
	return(out)
}

#' Wald tests for islptw objects
#' 
#' @param object An object of class \code{"islptw"}
#' @param est A vector indicating for which estimators to return a 
#' confidence interval. Possible estimators include the TMLE IPTW 
#' (\code{"islptw_tmle"}, recommended), the one-step IPTW 
#' (\code{"islptw_os"}, not recommended), the standard IPTW 
#' (\code{"iptw"}, recommended only for comparison to the other two estimators).
#' @param null The null hypothesis value(s). 
#' @param contrast This option specifies what parameter to return confidence intervals for.
#' If \code{contrast=NULL}, then test the null hypothesis that the covariate-adjusted
#' marginal means equal the value(s) specified in \code{null}. 
#' \code{contrast} can also be a numeric vector of ones, negative ones,
#' and zeros to define linear combinations of the various means (e.g., to estimate an average
#' treatment effect, see examples). In this case, we test the null hypothesis that the 
#' linear combination of means equals the value specified in \code{null}. 
#' \code{contrast} can also be a list with named functions
#' \code{f}, \code{h}, and \code{h_grad}. The function \code{f} takes
#' as input argument \code{eff} and specifies which transformation 
#' of the effect measure to test. The function \code{h} defines the contrast 
#' to be estimated and should take as input \code{est}, a vector
#' of the same length as \code{object$a_0}, and output the desired contrast. The function
#' \code{h_grad} is the gradient of the function \code{h(f())}. The function
#' computes a test of the null hypothesis that \code{h(f(object$est)) = null}. 
#' See examples. 
#' @param ... Other options (not currently used)
#' 
#' @importFrom stats pnorm
#' @export
#' @method wald_test islptw
#' @return An object of class \code{"ci.islptw"} with point estimates and
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
#' # get test that each mean = 0.5
#' test_mean <- wald_test(fit1, null = 0.5)
#' 
#' # get test that the ATE = 0 
#' ci_ATE <- ci(fit1, contrast = c(1,-1), null = 0)
#' 
#' # get test for risk ratio = 1 on log scale
#' myContrast <- list(f = function(eff){ log(eff) },
#'                    f_inv = function(eff){ exp(eff) }, # not necessary
#'                    h = function(est){ est[1]/est[2] },
#'                    h_grad =  function(est){ c(1/est[1],-1/est[2]) })
#' ci_RR <- ci(fit1, contrast = myContrast, null = 1)

wald_test.islptw <- function(object, est = c("islptw_tmle"), null = 0,
                             contrast = NULL, ...){
	if(any(est=="iptw")){
		stop("Theory does not support inference for naive IPTW with super learner.")
	}
	if(class(object) != "islptw"){
		stop("ci only works with islptw objects")
	}
	out <- vector(mode = "list", length = length(est))
	names(out) <- est
	# if no contrast then return an CI for each 
	# covariate-adjusted mean
	if(is.null(contrast)){
		# make sure null is input properly
		if(length(null) == 1){
			null <- rep(null, length(object$a_0))
		}else if(length(null) > length(object$a_0)){
			stop("length of null should be one or same length as object$a_0")
		}

		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = length(object$a_0), ncol = 2)
			for(j in seq_along(object$a_0)){
				zstat <- (object[[est[i]]]$est[j]-null[j])/sqrt(object[[est[i]]]$cov[j,j])
				pval <- 2*stats::pnorm(-abs(zstat))
				out[[i]][j,] <- c(zstat,pval)
			}
			row.names(out[[i]]) <- paste0(
              "H0: E[Y(",object$a_0,")]=",null
          	)
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}else if(is.numeric(contrast)){
		# check that contrast came in correctly
		if(length(contrast) != length(object$a_0)){
			stop("length of contrast vector not equal to length of a_0")
		}
		if(!all(contrast %in% c(-1,1,0))){
			stop("contrast should only be -1, 1, or 0")
		}
		if(length(null) > 1){
			stop("null length should be 1 if contrast is a linear combination")
		}
		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 2)
			g <- matrix(contrast, nrow = length(contrast))
			p <- matrix(object[[est[i]]]$est, nrow = length(contrast))
			v <- object[[est[i]]]$cov
			thisC <- t(g)%*%p
			thisSe <- sqrt(t(g)%*%v%*%g)
			zstat <- (thisC - null) / thisSe
			pval <- 2*stats::pnorm(-abs(zstat))
			out[[i]][1,] <- c(zstat, pval)
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
			thisName <- paste0("H0:",paste0(plusTerms,collapse = "+"),
			                   ifelse(length(indMinus)>0,"-",""),
			                   paste0(minusTerms,collapse="-"),"=",null)

			row.names(out[[i]]) <- thisName
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}else if(is.list(contrast)){
		if(!all(c("f","h","h_grad") %in% names(contrast))){
			stop("some function missing in contrast. see ?wald_test for help.")
		}
		if(length(null) > 1){
			stop("null length should be 1 if contrast is user-specified functions")
		}
		for(i in seq_along(est)){
			out[[i]] <- matrix(NA, nrow = 1, ncol = 2)
			thisC <- do.call(contrast$h,args=list(est = object[[est[i]]]$est))
			f_thisC <- do.call(contrast$f,args=list(eff = thisC))

			grad <- matrix(do.call(contrast$h_grad,
			                       args=list(est = object[[est[i]]]$est)
			               ),nrow = length(object$a_0))
			v <- object[[est[i]]]$cov
			thisSe <- sqrt(t(grad)%*%v%*%grad)
			zstat <- (f_thisC - null)/ thisSe
			pval <- 2*stats::pnorm(-abs(zstat))
			out[[i]][1,] <- c(zstat,pval)
			row.names(out[[i]]) <- paste0("H0: user contrast = ",null)
			colnames(out[[i]]) <- c("zstat","pval")
		}
	}
	class(out) <- "wald_test.islptw"
	return(out)
}