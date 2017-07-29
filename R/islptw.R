#' Compute asymptotically linear IPTW estimators with super learning 
#' for the propensity score
#' 
#' @param W A \code{data.frame} of named covariates
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes. 
#' @param a_0 A vector of treatment levels at which to compute the adjusted mean outcome. 
#' @param SL_g A vector of characters or a list describing the Super Learner library to be used 
#' for the propensity score. See \code{link{SuperLearner::SuperLearner}} for details.
#' @param SL_Qr A vector of characters or a list describing the Super Learner library to be used 
#' for the first reduced-dimension regression. 
#' @param glm_g A character describing a formula to be used in the call to \code{glm} for the propensity score. Ignored
#' if \code{SL_g!=NULL}.
#' @param glm_Qr A character describing a formula to be used in the call to \code{glm} for the first reduced-dimension regression. Ignored
#' if \code{SL_Qr!=NULL}. The formula should use the variable name \code{'gn'}.
#' @param maxIter A numeric that sets the maximum number of iterations the TMLE can perform in its fluctuation step.
#' @param tolIC A numeric that defines the stopping criteria based on the empirical mean
#' of the scores of the fluctuation submodels submodels. Setting to \code{"default"}
#' @param tolg A numeric indicating the minimum value for estimates of the propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param returnModels A boolean indicating whether to return model fits for the propensity score
#' and reduced-dimension regressions.
#' @return An object of class \code{"islptw"}.
#' \describe{
#'  \item{\code{islptw_tmle}}{A \code{list} of point estimates and 
#'        covariance matrix for the IPTW estimator based on a targeted propensity 
#' 		  score. }
#'  \item{\code{islptw_tmle_nuisance}}{A \code{list} of the final TMLE estimates of the
#'        propensity score (\code{$gnStar}) and reduced-dimension regression (\code{$QrnStar}) 
#' 		  evaluated at the observed data values.}
#'  \item{\code{islptw_os}}{A \code{list} of point estimates and covariance matrix for the
#' 		  one-step correct IPTW estimator.}
#'  \item{\code{islptw_os_nuisance}}{A \code{list} of the initial estimates of the
#'        propensity score and reduced-dimension regression evaluated at the observed data values.}
#'  \item{\code{iptw}}{A \code{list} of point estimates for the standard IPTW estimator. No
#' 		  estimate of the covariance matrix is provided because theory does not support asymptotic
#' 		  Normality of the IPTW estimator if super learning is used to estimate the propensity score.}
#'  \item{\code{aiptw}}{A \code{list} of doubly-robust point estimates and non-doubly-robust
#'        covariance matrix for the standard AIPTW estimator.}
#'  \item{\code{gcomp}}{A \code{list} of non-doubly-robust point estimates and non-doubly-robust
#'        covariance matrix for the standard G-computation estimator. If super learner is used
#'        there is no guarantee of correct inference for this estimator.}
#'  \item{\code{QnMod}}{The fitted object for the outcome regression. Returns \code{NULL}
#'        if \code{returnModels = FALSE}.}
#'  \item{\code{gnMod}}{The fitted object for the propensity score. Returns \code{NULL} if
#'        \code{returnModels = FALSE}.}
#'  \item{\code{QrnMod}}{The fitted object for the reduced-dimension regression that guards 
#'        against misspecification of the outcome regression. Returns \code{NULL} if 
#'        \code{returnModels = FALSE}.}
#'  \item{\code{a_0}}{The treatment levels that were requested for computation of 
#'        covariate-adjusted means.}
#' }
#' @importFrom plyr llply laply
#' @export
#' @examples 
#' # load super learner
#' library(SuperLearner)
#' # simulate data 
#' set.seed(123456)
#' n <- 200
#' W <- data.frame(W1 = runif(n), W2 = rnorm(n))
#' A <- rbinom(n,1,plogis(W$W1 - W$W2))
#' Y <- rbinom(n, 1, plogis(W$W1*W$W2*A))
#'
#' fit1 <- islptw(W = W, A = A, Y = Y, a_0 = c(1,0),
#'                SL_g=c("SL.glm","SL.mean","SL.step"),
#'                SL_Qr="SL.npreg")

islptw <- function(W, A, Y, a_0 = unique(A),
                   SL_g = NULL, 
                   glm_g = NULL,
                   SL_Qr = NULL,
                   glm_Qr = NULL,
                   returnModels = TRUE,
                   verbose = FALSE,
                   maxIter = 2, 
                   tolIC=1/length(Y), 
                   tolg=1e-2
                   ){

 	# estimate g
  	gnOut <- estimateG(A=A, W=W, tolg=tolg, verbose=verbose, 
  	                   returnModels=returnModels,SL_g=SL_g, 
  	                   glm_g=glm_g, a_0=a_0)
  	gn <- gnOut$est
 	

 	# compute iptw estimator
	psi.n <- mapply(a = split(a_0, 1:length(a_0)),g=gn, function(a,g){
		mean(as.numeric(A==a)/g * Y)
	})

    # estimate influence function
  	Dno <- mapply(a=split(a_0,1:length(a_0)),g=gn,psi=psi.n,FUN=function(a,g,psi){
    	as.numeric(A==a)/g * Y - psi
  	},SIMPLIFY=FALSE)

  	# estimate reduced dimension Q
  	# note that NULL is input to estimateQrn -- internally the function 
  	# assign Qn = 0 for all a_0 because estimateQrn estimates the regression
  	# of Y - Qn on gn (which is needed for drtmle), while here we just need
  	# the regression of Y on gn. 
    QrnOut <- estimateQrn(Y=Y, A=A, W=W, Qn=NULL, gn=gn, glm_Qr=glm_Qr, 
                          SL_Qr=SL_Qr, a_0=a_0,
                          returnModels = returnModels)
    Qrn <- QrnOut$est
    Dngo <- mapply(a=split(a_0,1:length(a_0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
      Qr/g * (as.numeric(A==a) - g)
    },SIMPLIFY=FALSE)
    PnDgn <- lapply(Dngo, mean)

    # one-step iptw estimator
    psi.o <- mapply(a=psi.n,b=PnDgn,SIMPLIFY=FALSE,
    	            FUN=function(a,b){a-b}) 

    # targeted g estimator
    gnStar <- gn
    QrnStar <- Qrn
  	PnDgnStar <- Inf
    eps <- list(Inf)
    ct <- 0
    # fluctuate
    while(max(abs(unlist(PnDgnStar))) > tolIC & ct < maxIter){
      ct <- ct + 1
    
      # fluctuate gnStar
      gnStarOut <- fluctuateG(Y=Y, A=A, W=W, a_0=a_0, tolg=tolg, 
                              gn=gnStar, Qrn=QrnStar)
      gnStar <- plyr::llply(gnStarOut, function(x){unlist(x$est)})
      eps <- plyr::laply(gnStarOut, function(x){x$eps})
      # re-estimate reduced dimension regression
      QrnStarOut <- estimateQrn(Y=Y, A=A, W=W, a_0=a_0, Qn=NULL, gn=gnStar, 
                                glm_Qr=glm_Qr, SL_Qr=SL_Qr, 
                                returnModels = returnModels)
      QrnStar <- QrnStarOut$est
      # compute influence function for fluctuated estimators
      DngoStar <- mapply(a=split(a_0,1:length(a_0)),Qr=QrnStar,g=gnStar,FUN=function(a,Qr,g){
        Qr/g * (as.numeric(A==a) - g)
      }, SIMPLIFY=FALSE)
      PnDgnStar <- lapply(DngoStar, mean)
      if(verbose){
        cat("TMLE Iteration", ct, "=", round(unlist(eps),5), "\n")
        cat("Mean of IC       =", round(unlist(PnDgnStar), 10),"\n")
      }
    }

    # compute final tmle-iptw estimate
 	# compute iptw estimator
	psi.nStar <- mapply(a = split(a_0, 1:length(a_0)),g=gnStar, function(a,g){
		mean(as.numeric(A==a)/g * Y)
	})

	# compute variance estimators 
	# original influence function
  	DnoStar <- mapply(a=split(a_0,1:length(a_0)),g=gnStar,
  	                  psi=psi.nStar,FUN=function(a,g,psi){
  	                  	as.numeric(A==a)/g * Y - psi
  	                  },SIMPLIFY=FALSE)
  	# new influence function in DngoStar already
    DnoStarMat <- matrix(unlist(DnoStar) - unlist(DngoStar), 
                         ncol=length(Y), nrow=length(a_0), byrow = TRUE)
    cov.t <- (DnoStarMat - mean(DnoStarMat))%*%t((DnoStarMat - mean(DnoStarMat)))/(length(Y)^2)
    DnoMat <- matrix(unlist(Dno) - unlist(DngoStar), 
                     ncol = length(Y), nrow = length(a_0), byrow = TRUE)
    cov.os <- (DnoMat - mean(DnoMat))%*%t((DnoMat - mean(DnoMat)))/(length(Y)^2)

    # output
    out <- list(islptw_tmle = list(est = unlist(psi.nStar), cov = cov.t),
                islptw_tmle_nuisance = list(gn = gnStar, QrnStar = QrnStar),
                islptw_os = list(est = unlist(psi.o), cov = cov.os),
                islptw_os_nuisance = list(gn = gn, Qrn = Qrn),
                iptw = list(est = unlist(psi.n)),
                gnMod = NULL, QrnMod = NULL, a_0 = a_0)
  	if(returnModels){
  	 out$gnMod <- gnOut$fm
   	 out$QrnMod <- QrnStarOut$fm
  	}
  	class(out) <- "islptw"
  	return(out)
}