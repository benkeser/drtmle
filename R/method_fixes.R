#' Modification of convex combination least squares method 
#' for \code{SuperLearner} that doesn't fail with duplicated predictors.
#' Issue reported to \code{SuperLearner} maintainers. This function will
#' be deprecated when a more robust fix in the \code{SuperLearner} package
#' is implemented.
#' @export
method.CC_LS_mod <- function()
{
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, 
        ...) {
        cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x - 
            Y)^2))
        names(cvRisk) <- libraryNames
        compute <- function(x, y, wt = rep(1, length(y))) {
            wX <- sqrt(wt) * x
            wY <- sqrt(wt) * y
            D <- crossprod(wX)
            d <- crossprod(wX, wY)
            A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
            bvec <- c(1, rep(0, ncol(wX)))
            fit <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, 
                bvec = bvec, meq = 1)
            invisible(fit)
        }
        colDup <- which(duplicated(round(Z, 8), MARGIN = 2))
        modZ <- Z
        if(length(colDup) > 0){
        	warning(paste0("Algorithm ", colDup, " is duplicated. Setting weight to 0."))
        	modZ <- modZ[,-colDup]
        }
        fit <- compute(x = modZ, y = Y, wt = obsWeights)
        coef <- fit$solution
        if(length(colDup) > 0){
	        ind <- c(seq_along(coef),colDup-0.5)
			coef <- c(coef,rep(0,length(colDup)))
	        coef <- coef[order(ind)]
		}
        if (any(is.na(coef))) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] = 0
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        if (!sum(coef) > 0) 
            warning("All algorithms have zero weight", call. = FALSE)
        list(cvRisk = cvRisk, coef = coef, optimizer = fit)
    }
    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }
    out <- list(require = "quadprog", computeCoef = computeCoef, 
        computePred = computePred)
    invisible(out)
}

#' Modification of convex combination negative log likelihood method
#' for \code{SuperLearner} that doesn't fail with duplicated predictors.
#' Issue reported to \code{SuperLearner} maintainers. This function will
#' be deprecated when a more robust fix in the \code{SuperLearner} package
#' is implemented.
#' @importFrom SuperLearner trimLogit
#' @importFrom stats plogis
#' @export
method.CC_nloglik_mod <- function () 
{
    computePred = function(predY, coef, control, ...) {
        if (sum(coef != 0) == 0) {
            stop("All metalearner coefficients are zero, cannot compute prediction.")
        }
        stats::plogis(SuperLearner::trimLogit(predY[, coef != 0], trim = control$trimLogit) %*% 
            matrix(coef[coef != 0]))
    }
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, 
        verbose, ...) {
        colDup <- which(duplicated(round(Z, 8), MARGIN = 2))
        modZ <- Z
        if(length(colDup) > 0){
            warning(paste0("Algorithm ", colDup, " is duplicated. Setting weight to 0."))
            modZ <- modZ[,-colDup]
        }
        modlogitZ = SuperLearner::trimLogit(modZ, control$trimLogit)
        logitZ = SuperLearner::trimLogit(Z, control$trimLogit)
        cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights * 
            ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x, log.p = TRUE, 
                lower.tail = FALSE))))
        names(cvRisk) <- libraryNames
        obj_and_grad <- function(y, x, w = NULL) {
            y <- y
            x <- x
            function(beta) {
                xB <- x %*% cbind(beta)
                loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 - 
                  y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
                if (!is.null(w)) 
                  loglik <- loglik * w
                obj <- -2 * sum(loglik)
                p <- stats::plogis(xB)
                grad <- if (is.null(w)) 
                  2 * crossprod(x, cbind(p - y))
                else 2 * crossprod(x, w * cbind(p - y))
                list(objective = obj, gradient = grad)
            }
        }

        lower_bounds = rep(0, ncol(modZ))
        upper_bounds = rep(1, ncol(modZ))
        if (anyNA(cvRisk)) {
            upper_bounds[is.na(cvRisk)] = 0
        }
        r <- nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)), eval_f = obj_and_grad(Y, 
            modlogitZ), lb = lower_bounds, ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) - 
            1), eval_jac_g_eq = function(beta) rep(1, length(beta)), 
            opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
        if (r$status < 1 || r$status > 4) {
            warning(r$message)
        }
        coef <- r$solution
        if(length(colDup) > 0){
            ind <- c(seq_along(coef),colDup - 0.5)
            coef <- c(coef,rep(0,length(colDup)))
            coef <- coef[order(ind)]
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
        return(out)
    }
    list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}