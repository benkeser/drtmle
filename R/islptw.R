#' Compute asymptotically linear IPTW estimators with super learning 
#' for the propensity score
#' 
#' @param W A \code{data.frame} of named covariates
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes.
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed) 
#' @param a_0 A vector of treatment levels at which to compute the adjusted mean outcome. 
#' @param stratify A \code{boolean} indicating whether to estimate the missing outcome regression separately
#' for observations with different levels of \code{A} (if \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
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
#' @param cvFolds A numeric equal to the number of folds to be used in cross-validated fitting of 
#' nuisance parameters. If \code{cvFolds = 1}, no cross-validation is used.
#' @param ... Other options (not currently used).
#' @param parallel A boolean indicating whether to use \code{foreach}
#' to estimate nuisance parameters in parallel. Only useful if there is a registered parallel
#' backend (see examples) and \code{cvFolds > 1}.
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

islptw <- function(W, A, Y, 
                   DeltaY = as.numeric(!is.na(Y)), 
                   DeltaA = as.numeric(!is.na(A)), 
                   stratify = FALSE, 
                   a_0 = unique(A[!is.na(A)]),
                   SL_g = NULL, 
                   glm_g = NULL,
                   SL_Qr = NULL,
                   glm_Qr = NULL,
                   returnModels = TRUE,
                   verbose = FALSE,
                   maxIter = 2, 
                   tolIC=1/length(Y), 
                   tolg=1e-2,
                   cvFolds = 1, parallel = FALSE,
                   ... 
                   ){
  # if cvFolds non-null split data into cvFolds pieces
  n <- length(Y)
  if(cvFolds != 1){
    validRows <- split(sample(1:n), rep(1:cvFolds, length = n))       
    ordVR <- order(unlist(validRows))
  }else{
    validRows <- list(NULL)
    ordVR <- 1:n
  }

  #-------------------------------
  # estimate propensity score
  #-------------------------------
  if(!parallel){
    gnOut <- lapply(X = validRows, FUN = estimateG,
                    A=A, W=W, DeltaA = DeltaA, DeltaY = DeltaY, 
                    tolg=tolg, verbose=verbose, 
                    returnModels=returnModels,SL_g=SL_g,
                    glm_g=glm_g, a_0=a_0,stratify = stratify)
  }else{
    gnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
      estimateG(A = A, W = W, DeltaA = DeltaA, DeltaY = DeltaY, 
                tolg = tolg, verbose = verbose,stratify = stratify,
                returnModels = returnModels, SL_g = SL_g,
                glm_g = glm_g, a_0 = a_0, validRows = validRows[[v]])
    }
  }
  # re-order predictions
  gnValid <- unlist(gnOut, recursive = FALSE, use.names = FALSE)
  gnUnOrd <- do.call(Map, c(c, gnValid[seq(1, length(gnValid), 2)]))
  gn <- lapply(gnUnOrd, function(x){ x[ordVR] })
  # obtain list of propensity score fits
  gnMod <- gnValid[seq(2, length(gnValid), 2)]
  # TO DO: Add reasonable names to gnMod?

 	# compute iptw estimator
	psi.n <- mapply(a = split(a_0, 1:length(a_0)),g=gn, function(a,g){
    modA <- A; modA[is.na(A)] <- -999
    modY <- Y; modY[is.na(Y)] <- -999
		mean(as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * modY)
	})

  # estimate influence function
	Dno <- mapply(a=split(a_0,1:length(a_0)),g=gn,psi=psi.n,FUN=function(a,g,psi){
    modA <- A; modA[is.na(A)] <- -999
    modY <- Y; modY[is.na(Y)] <- -999
  	as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * modY - psi
	},SIMPLIFY=FALSE)

  #-------------------------------------
	# estimate reduced dimension Q
  #-------------------------------------
	# note that NULL is input to estimateQrn -- internally the function 
	# assign Qn = 0 for all a_0 because estimateQrn estimates the regression
	# of Y - Qn on gn (which is needed for drtmle), while here we just need
	# the regression of Y on gn. 
  if(!parallel){
    QrnOut <- lapply(X = validRows, FUN = estimateQrn, 
                   Y=Y, A=A, W=W, DeltaA = DeltaA, DeltaY = DeltaY, 
                   Qn=NULL, gn=gn, glm_Qr=glm_Qr, 
                   SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels)  
  }else{
    QrnOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
      estimateQrn(Y=Y, A=A, W=W, DeltaA = DeltaA, DeltaY = DeltaY, 
                  Qn=NULL, gn=gn, glm_Qr=glm_Qr, 
                  SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels,
                  validRows = validRows[[v]])
    }
  }

  # re-order predictions
  QrnValid <- unlist(QrnOut, recursive = FALSE, use.names = FALSE)
  QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid),2)]))
  Qrn <- lapply(QrnUnOrd, function(x){ x[ordVR] })
  # obtain list of propensity score fits
  QrnMod <- QrnValid[seq(2,length(QrnValid),2)]
  # TO DO: Add reasonable names to QrnMod?

  Dngo <- mapply(a=split(a_0,1:length(a_0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
    modA <- A; modA[is.na(A)] <- -999
    Qr/g * (as.numeric(modA == a & DeltaA == 1 & DeltaY == 1) - g)
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
    gnStarOut <- fluctuateG(Y=Y, A=A, W=W, 
                            DeltaA = DeltaA, DeltaY = DeltaY, 
                            a_0=a_0, tolg=tolg, 
                            gn=gnStar, Qrn=QrnStar)
    gnStar <- plyr::llply(gnStarOut, function(x){unlist(x$est)})
    eps <- plyr::laply(gnStarOut, function(x){x$eps})
    # re-estimate reduced dimension regression
    if(!parallel){
      QrnStarOut <- lapply(X = validRows, FUN = estimateQrn, 
                   Y=Y, A=A, W=W, DeltaA = DeltaA, DeltaY = DeltaY, 
                   Qn=NULL, gn=gnStar, glm_Qr=glm_Qr, 
                   SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels)
    }else{
      QrnStarOut <- foreach::foreach(v = 1:cvFolds, .packages = "SuperLearner") %dopar% {
        estimateQrn(Y=Y, A=A, W=W, DeltaA = DeltaA, DeltaY = DeltaY, 
                    Qn=NULL, gn=gnStar, glm_Qr=glm_Qr, 
                    SL_Qr=SL_Qr, a_0=a_0,returnModels = returnModels,
                    validRows = validRows[[v]])
      }
    }
    # re-order predictions
    QrnValid <- unlist(QrnStarOut, recursive = FALSE, use.names = FALSE)
    QrnUnOrd <- do.call(Map, c(c, QrnValid[seq(1,length(QrnValid),2)]))
    QrnStar <- lapply(QrnUnOrd, function(x){ x[ordVR] })
    # obtain list of propensity score fits
    QrnMod <- QrnValid[seq(2,length(QrnValid),2)]

    # compute influence function for fluctuated estimators
    DngoStar <- mapply(a=split(a_0,1:length(a_0)),
                       Qr=QrnStar,g=gnStar,
                       FUN=function(a,Qr,g){
      modA <- A; modA[is.na(A)] <- -999
      Qr/g * (as.numeric(modA == a & DeltaA == 1 & DeltaY == 1) - g)
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
    modA <- A; modA[is.na(A)] <- -999
    modY <- Y; modY[is.na(Y)] <- -999
  	mean(as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * modY)
  })

  # compute variance estimators 
  # original influence function
	DnoStar <- mapply(a=split(a_0,1:length(a_0)),g=gnStar,
	                  psi=psi.nStar,FUN=function(a,g,psi){
                      modA <- A; modA[is.na(A)] <- -999
                      modY <- Y; modY[is.na(Y)] <- -999
	                  	as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * modY - psi
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
	 out$gnMod <- gnMod
 	 out$QrnMod <- QrnMod
	}
	class(out) <- "islptw"
	return(out)
}