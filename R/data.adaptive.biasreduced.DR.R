#' data.adaptive.biasreduced.DR
#' Vermeulen estimator from IJB paper
#' @param R Missingness indicator
#' @param Y Outcome
#' @param cov Covariates
#' @param type.initQ How to estimate outcome regression (\code{"par","SL","npreg"})
#' @param zeta Truncation level for logit transform of Qn
#' @param fluc How to fit fluctuation submodel 
#' @param alpha Truncation level for propensity score gn
#' @param psi.tilde Null hypothesis value
#' 
#' @importFrom stats qnorm coef pnorm sd
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
                SL.library <- c("SL.gbm.caret1","SL.step.interaction","SL.glm")
                initQ <- (SuperLearner(Y=Ym.star,X=dat.cov.m,
                                      newX=dat.cov,verbose = FALSE,
                                      SL.library=SL.library,
                                      method="method.NNLS")$SL.predict -a)/(b-a)
            }else if(type.initQ=="npreg"){
                dat.cov.m <- dat.cov[R==1,,drop=FALSE]
                fm <- do.call("SL.npreg",args=list(Y=Y[R==1], X=dat.cov.m,obsWeights=rep(1,length(Y[R==1])),
                                                 newX=dat.cov, family=data.frame(family="gaussian")))
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
                                   family="binomial",offset=logit(initQ.trunc),
                                   subset=(R==1))
                flucQ <- expit(logit(initQ.trunc)+
                                   as.vector(stats::coef(fluctuationQ)%*%t(w.cov)))
            }
            else if(fluc=="weighted"){
                fluctuationQ <-glm(Y.star~cov,
                                   family="binomial",offset=logit(initQ.trunc),
                                   subset=(R==1), weights=(1-ps.par)/ps.par)
                flucQ <- expit(logit(initQ.trunc)+
                   as.vector(stats::coef(fluctuationQ)%*%t(int.cov)))
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
        se.psi <- (b-a)*stats::sd(U(R=R,Y=Y.star,
                             outcome=flucQ,ps=ps.par))/sqrt(n)
        # 95% confidence interval
        ci.psi <- psi+c(-1,1)*stats::qnorm(1-alpha/2)*se.psi
        # Wald test statistic
        W <- (psi-psi.tilde)/se.psi
        # p-value Wald test
        p.value <- 2*stats::pnorm(abs(W),lower.tail=FALSE)
        return(list(est=psi,se=se.psi,ci=ci.psi,
                    Wald.statistic=W,p.value=p.value))
    }



