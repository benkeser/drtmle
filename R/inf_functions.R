#' Evaluate usual efficient influence function
#' 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes. 
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qn List of estimated outcome regression evaluated at observations
#' @param gn List of estimated propensity scores evaluated at observations
#' @param psi_n List of estimated ATEs
#' @param a_0 Vector of values to return marginal mean
#' 

eval_Dstar <- function(A, Y, DeltaY, DeltaA, Qn, gn, psi_n, a_0){
	return(mapply(a=split(a_0,1:length(a_0)),Q=Qn,g=gn,p=psi_n,FUN=function(a,Q,g,p){
		# in order for influence function computations to compute properly
  		# replace missing A and Y by arbitrary numerics.
  		# Note that when these are missing, they are always getting multiplied
  		# by 0, but R return 0*NA = NA for some reason and this is a hacky fix
  		# to get around that. 
    	modA <- A; modA[is.na(A)] <- -999
    	modY <- Y; modY[is.na(Y)] <- -999
    	as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * (modY - Q) + Q - p
  	},SIMPLIFY=FALSE))
}

#' Evaluate extra piece of efficient influence function resulting from
#' misspecification of outcome regression
#' 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qrn List of estimated reduced-dimension outcome regression evaluated at observations
#' @param gn List of estimated propensity scores evaluated at observations
#' @param a_0 Vector of values to return marginal mean
#' 

eval_Dstar_g <- function(A, DeltaY, DeltaA, Qrn, gn, a_0){
	return(mapply(a=split(a_0,1:length(a_0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
      modA <- A; modA[is.na(A)] <- -999
      Qr/g * (as.numeric(modA == a & DeltaA == 1 & DeltaY == 1) - g)
    },SIMPLIFY=FALSE))
}


#' Evaluate extra piece of efficient influence function resulting from
#' misspecification of propensity score
#' 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes. 
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qn List of estimated outcome regression evaluated at observations
#' @param gn List of estimated propensity scores evaluated at observations
#' @param grn List of estimated reduced-dimension propensity scores evaluated at observations
#' @param a_0 Vector of values to return marginal mean
#' @param reduction A character equal to \code{"univariate"} for a univariate misspecification correction or \code{"bivariate"}
#' for the bivariate version.
#' 
eval_Dstar_Q <- function(A, Y, DeltaY, DeltaA, Qn, gn, grn, a_0, reduction){
	if(reduction=="univariate"){
      return(mapply(a=split(a_0,1:length(a_0)),Q=Qn,gr=grn,FUN=function(a,Q,gr){
        modA <- A; modA[is.na(A)] <- -999
        modY <- Y; modY[is.na(Y)] <- -999
        as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/gr$grn2 * gr$grn1 * (modY-Q)
      },SIMPLIFY=FALSE))
    }else if(reduction=="bivariate"){
      return(mapply(a=split(a_0,1:length(a_0)),Q=Qn,g=gn, gr=grn,FUN=function(a,Q,gr,g){
        modA <- A; modA[is.na(A)] <- -999
        modY <- Y; modY[is.na(Y)] <- -999
        as.numeric(modA==a & DeltaA == 1 & DeltaY == 1)/gr$grn2 * (gr$grn2 - g)/g * (modY-Q)
      },SIMPLIFY=FALSE))
    }
}

#' Evaluate usual influence function of IPTW
#' 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param Y A numeric of continuous or binary outcomes. 
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param gn List of estimated propensity scores evaluated at observations
#' @param a_0 Vector of values to return marginal mean
#' @param psi_n List of estimated ATEs

eval_Diptw <- function(A, Y, DeltaA, DeltaY, gn, psi_n, a_0){
	return(mapply(a=split(a_0,1:length(a_0)),g=gn,psi=psi_n,FUN=function(a,g,psi){
    	modA <- A; modA[is.na(A)] <- -999
    	modY <- Y; modY[is.na(Y)] <- -999
  		as.numeric(modA == a & DeltaA == 1 & DeltaY == 1)/g * modY - psi
	},SIMPLIFY=FALSE))
}


#' Evaluate extra piece of the influence function for the IPTW
#' 
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or 1)
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if missing 1 if observed)
#' @param Qrn List of estimated reduced-dimension outcome regression evaluated at observations
#' @param gn List of estimated propensity scores evaluated at observations
#' @param a_0 Vector of values to return marginal mean
#' 
eval_Diptw_g <- function(A, DeltaA, DeltaY, Qrn, gn, a_0){
	return(mapply(a=split(a_0,1:length(a_0)),Qr=Qrn,g=gn,FUN=function(a,Qr,g){
	    modA <- A; modA[is.na(A)] <- -999
	    Qr/g * (as.numeric(modA == a & DeltaA == 1 & DeltaY == 1) - g)
  	},SIMPLIFY=FALSE))
}