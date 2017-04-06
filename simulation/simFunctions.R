#' makeData
#' 
#' Make Kang and Schafer data set of size n

makeData <- function(n){
    Z1 <- rnorm(n)
    Z2 <- rnorm(n)
    Z3 <- rnorm(n)
    Z4 <- rnorm(n)
    
    W <- data.frame(
        W1 = exp(Z1/2),
        W2 = Z2 / (1+exp(Z1)) + 10,
        W3 = (Z1*Z3/25 + 0.6)^3,
        W4 = (Z2 + Z4 + 20)^2
    )
    
    Q0 <- 210 + 27.4*Z1 + 13.7*(Z2 + Z3+ Z4) 
    Y <- rnorm(n, Q0, 1)
    
    g0 <- plogis(-Z1 + 0.5*Z2 - 0.25*Z3 -0.1*Z4)
    A <- rbinom(n, 1, g0)
    
    return(list(W=W, A=A, Y=Y))
}

do.one <- function(n, seed, beta, guard){
  set.seed(seed)
    n <- 500
    library(plyr)
    library(np)
    library(SuperLearner)
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  Z3 <- rnorm(n)
  Z4 <- rnorm(n)
  # Z5 <- rnorm(n)
  # Z6 <- rnorm(n)
  # Z7 <- rnorm(n)
  # Z8 <- rnorm(n)
  Z <- eval(parse(text=paste0("data.frame(",paste0("Z",1:8,collapse=","),")")))
  
  W <- data.frame(
    W1 = exp(Z1/2),
    W2 = Z2 / (1+exp(Z1)) + 10,
    W3 = (Z1*Z3/25 + 0.6)^3,
    W4 = (Z2 + Z4 + 20)^2
  )
  #   W5 = log(abs(Z5 + 5) + 0.5),
  #   W6 = Z6 / (1+exp(Z7)) + 1,
  #   W7 = exp(Z8)/(1+exp(Z7)),
  #   W8 = (Z8 - Z7)/(Z6+4)
  # )

  Q0 <- 210 + 27.4*Z1 + 13.7*(Z2 + Z3+ Z4) #+ Z5+ Z6+ Z7+ Z8) 
  Y <- rnorm(n, Q0, 1)

  g0 <- plogis(-Z1 + 0.5*Z2 - 0.25*Z3 -0.1*Z4)
  A <- rbinom(n, 1, g0)

  set.seed(1234)
  out <- drtmle(Y=Y, A=A, W=W, family=gaussian(),
                   a0=1,
                   guard=c("Q","g"),
                   libraryQ=c("SL.hal"),
                   libraryg=c("SL.hal"),
                   librarygr=c("SL.npreg"),
                   libraryQr=c("SL.npreg"),
                   reduction="univariate",
                   glmQ=NULL,
                   glmg=NULL,
                   glmgr=NULL,
                   glmQr=NULL, tolEps=1e-4, tolIC="default", 
                   tolg=0.025
  )

  
  # new vermeulen
  out2 <- data.adaptive.biasreduced.DR(
      R=A,Y=Y,cov=data.matrix(W),type.initQ=c("SL"),
      zeta=0.005,fluc=c("unweighted"),
      alpha=0.05,psi.tilde=0)
  
  # cao
  out3 <- cao.dr(R=A,Y=Y,cov=W,nBoot=100)
  
  # original vermeulen
  out4 <- m.biasreducedDR.identity(R=A,Y=Y,cov=data.matrix(W))
  
  # tan
  library(iWeigReg)
  Qmod <- glm(Y[A==1] ~ ., data = W[A==1,])
  Qn <- predict(Qmod, type="response", newdata=W)
  gmod <- glm(A ~ ., data=W, family="binomial")
  gn <- predict(gmod, type="response")
  
  out5 <- mn.clik(y = Y, tr = A, p = gn, g = cbind(1, Qn), X = model.matrix(gmod))
  
          
  return(
    list(seed=seed, n=n, beta=beta, out=out)
  )

}