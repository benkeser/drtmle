#' data.adaptive.biasreduced.DR
#' 
#' Vermeulen estimator from IJB paper
#' @export 

data.adaptive.biasreduced.DR <-
    function(R,Y,cov,type.initQ=c("par","SL"),
             zeta=0.005,fluc=c("unweighted","weighted"),
             alpha=0.05,psi.tilde=0){
        
        expit <- function(x){1/(1+exp(-x))}
        logit <- function(x){log(x/(1-x))}
        n <- length(R)
        dat.cov <- data.frame(cov)
        int.cov <- cbind(rep(1,n),cov)
        colnames(dat.cov)<-
            paste("cov.",1:dim(cov) [2],sep="")
        # propensity score
        mler <- glm(R~cov,family="binomial")
        ps.par <- predict(mler,type="response")
        # initial conditional mean outcome
        a <- min(Y[R==1]) - 0.1*abs(min(Y[R==1]))
        b <- max(Y[R==1]) + 0.1*abs(max(Y[R==1]))
        Y.star <- (Y-a)/(b-a)
        {
            if(type.initQ=="par") {
                mley <- lm(Y.star~cov,subset=(R==1))
                initQ <- predict(mley,newdata=dat.cov)
            }
            else if(type.initQ=="SL"){
                Ym.star <- Y.star[R==1]
                dat.cov.m <- dat.cov[R==1,,drop=FALSE]
                # SL.library <- c("SL.glm", "SL.randomForest",
                #                 "SL.gam","SL.polymars", "SL.mean")
                # initQ <- SuperLearner(Y=Ym.star,X=dat.cov.m,
                #                       newX=dat.cov,verbose = FALSE,
                #                       SL.library=SL.library,
                #                       method="method.NNLS")$SL.predict
                fm <- do.call("SL.hal",args=list(Y=Y[R==1], X=dat.cov.m,
                                                 newX=dat.cov, family=gaussian()))
                initQ <- (fm$pred-a)/(b-a)
            }
        }
        initQ.trunc <- ifelse(initQ < zeta,zeta,
                              ifelse(initQ > 1-zeta,1-zeta,initQ))
        # fluctuation
        {
            if(fluc=="unweighted"){
                w.cov <- (1-ps.par)/ps.par*int.cov
                fluctuationQ <- glm(Y.star ~ -1 + w.cov,
                                   family=binomial,offset=logit(initQ.trunc),
                                   subset=(R==1))
                flucQ <- expit(logit(initQ.trunc)+
                                   as.vector(coef(fluctuationQ)%*%t(w.cov)))
            }
            else if(fluc=="weighted"){
                fluctuationQ <-glm(Y.star~cov,
                                   family=binomial,offset=logit(initQ.trunc),
                                   subset=(R==1), weights=(1-ps.par)/ps.par)
                flucQ <- expit(logit(initQ.trunc)+
                   as.vector(coef(fluctuationQ)%*%t(int.cov)))
            }
        }
        # doubly robust estimator
        U <- function(R,Y,outcome,ps){
            outcome+R/ps*(Y-outcome)
        }
        est.trunc <- mean(U(R=R,Y=Y.star,
                            outcome=flucQ,ps=ps.par))
        psi <- (b-a)*est.trunc+a
        # standard error
        se.psi <- (b-a)*sd(U(R=R,Y=Y.star,
                             outcome=flucQ,ps=ps.par))/sqrt(n)
        # 95% confidence interval
        ci.psi <- psi+c(-1,1)*qnorm(1-alpha/2)*se.psi
        # Wald test statistic
        W <- (psi-psi.tilde)/se.psi
        # p-value Wald test
        p.value <- 2*pnorm(abs(W),lower.tail=FALSE)
        return(list(est=psi,se=se.psi,ci=ci.psi,
                    Wald.statistic=W,p.value=p.value))
    }


#' m.biasreducedDR.identity
#' 
#' Code to implement Vermeulen non-data-adaptive estimator
#' @export


m.biasreducedDR.identity<-function(R,Y,cov){
  n<-length(R)
  int.cov<-cbind(rep(1,n),cov)
  expit<-function(x) exp(x)/(1+exp(x))
  U <- function(R,Y,X,gamma,beta){
    (R/expit(gamma%*%t(X))*(Y-beta%*%t(X))+beta%*%t(X))
  }
  min.Uint<-function(gamma){
    -mean((-R*exp(-gamma%*%t(int.cov)))+(-(1-R)*(gamma%*%t(int.cov))))
  }
  init.gamma<-coef(glm(R~cov,family="binomial"))
  sol<-nlm(min.Uint,init.gamma)
  gamma.BR<-sol$estimate
  weight<-as.vector(1/exp(gamma.BR%*%t(int.cov)))
  beta.BR<-coef(lm(Y ~ -1+int.cov,subset=(R==1),weights=weight))
  mn.Y<-mean(U(R,Y,int.cov,gamma.BR,beta.BR))
  se.mn.Y<-sd(U(R,Y,int.cov,gamma.BR,beta.BR))/sqrt(n)
  return(list(mn.Y=mn.Y,se.mn.Y=se.mn.Y,gamma.BR=gamma.BR,beta.BR=beta.BR))
}


#' getCaoEst
#' 
#' Code to implement Cao 2009 estimator
getCaoEst <- function(R,Y,cov,family){
    # fit outcome regression
    Qmod <- glm(paste0("Y ~", paste0(colnames(cov),collapse="+")),
        data = data.frame(Y, cov)[R==1,],
        family = family)

    Qn <- predict(Qmod, newdata=data.frame(Y,cov), type="response")

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

    tmp <- optim(rep(0,ncol(cov)+2), fn = negLogLik, R = R, cov = cov,
     control=list(maxit=1e2))
    p <- tmp$par
    delta <- p[1]; gamma <- matrix(p[2:length(p)],ncol=1)
    X <- data.matrix(cbind(rep(1,length(R)), cov))
    gn <- 1-exp(delta + X%*%gamma)/(1 + exp(X%*%gamma))
    if(any(gn > 1) | any(gn < 0)){
        suppressWarnings(
        tmp2 <- constrOptim.nl(tmp$par, fn = negLogLik, hin = constraint, R = R, cov = cov)
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

#' cao.dr 
#' 
#' Compute the Cao 2009 estimator
#' @export 
cao.dr <- function(R,Y,cov, family = "gaussian",nBoot=500){
    # get estimate
    est <- getCaoEst(R=R,Y=Y,cov=cov,family=family)
    # get bootstrap CI
    doOneBoot <- function(){
        samp <- sample(1:length(R), replace = TRUE)
        tmpEst <- getCaoEst(R=R[samp],Y=Y[samp],cov=cov[samp,],
            family=family)
        return(tmpEst)
    }
    estVec <- replicate(nBoot, doOneBoot())
    ci <- quantile(estVec, p = c(0.025, 0.975))
    return(list(est = est, ci = ci))
}



