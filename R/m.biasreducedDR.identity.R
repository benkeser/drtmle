#' m.biasreducedDR.identity
#' 
#' @param R Missingness indicator
#' @param Y Outcome
#' @param cov Covariates
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

