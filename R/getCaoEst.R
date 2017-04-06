
#' getCaoEst
#' 
#' Code to implement Cao 2009 estimator
#' @param R Missing indicator
#' @param Y outcome
#' @param cov Covariates
#' @param family Character for family of outcome
#' 
#' @importFrom stats glm predict optim
#' @importFrom alabama constrOptim.nl
#' 
getCaoEst <- function(R,Y,cov,family){
    # fit outcome regression
    Qmod <- stats::glm(paste0("Y ~", paste0(colnames(cov),collapse="+")),
        data = data.frame(Y, cov)[R==1,],
        family = family)

    Qn <- stats::predict(Qmod, newdata=data.frame(Y,cov), type="response")

    # fit "enhanced propensity" regression
    negLogLik <- function(pars,R,cov){
    delta <- pars[1]; gamma <- matrix(pars[2:length(pars)],ncol=1)
    X <- data.matrix(cbind(rep(1,length(R)), cov))
    g <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    return(-sum(R*log(g) + (1-R)*log(1-g)))
    }

    constraint <- function(pars, R,cov){
     delta <- pars[1]; gamma <- matrix(pars[2:length(pars)],ncol=1)
     X <- data.matrix(cbind(rep(1,length(R)), cov))
     g <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
     # first that all probs are between 0 and 1
     c1 <- as.numeric(all(g < 1) & all(g > 0))
     # now that sum of inverse weights sum to n
     c2 <- as.numeric(sum(R/g) == length(R))    
     return(c1)
    }

    tmp <- stats::optim(rep(0,ncol(cov)+2), fn = negLogLik, R = R, cov = cov,
     control=list(maxit=1e2))
    p <- tmp$par
    delta <- p[1]; gamma <- matrix(p[2:length(p)],ncol=1)
    X <- data.matrix(cbind(rep(1,length(R)), cov))
    gn <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    if(any(gn > 1) | any(gn < 0)){
        suppressWarnings(
        tmp2 <- alabama::constrOptim.nl(tmp$par, fn = negLogLik, hin = constraint, R = R, cov = cov)
        )
        p <- tmp2$par
        delta <- p[1]; gamma <- matrix(p[2:length(p)],ncol=1)
        X <- data.matrix(cbind(rep(1,length(R)), cov))
        g2 <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    }
    
    est <- mean(
        R*Y/gn - (R - gn)/gn * Qn
    )
    return(est)
}
