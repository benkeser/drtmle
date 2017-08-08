#--------------------------------------------------
# Print methods
#--------------------------------------------------
#' Print the output of a \code{"drtmle"} object.
#' 
#' @param x A \code{"drtmle"} object
#' @param ... Other arguments (not used)
#' @export
#' @method print drtmle

print.drtmle <- function(x,...){
  tmp <- list(est = cbind(x$drtmle$est),
  	          cov = x$drtmle$cov)
  row.names(tmp$est) <- x$a_0
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if(length(x$a_0) <= 4){
  	print(tmp)
  }else{
  	tmp$cov <- diag(tmp$cov)
  }
}

#' Print the output of a \code{"islptw"} object.
#' 
#' @param x A \code{"islptw"} object.
#' @param ... Other arguments (not used)
#' @export
#' @method print islptw

print.islptw <- function(x,...){
  tmp <- list(est = cbind(x$islptw_tmle$est),
  	          cov = x$islptw_tmle$cov)
  row.names(tmp$est) <- x$a_0
  row.names(tmp$cov) <- colnames(tmp$cov) <- x$a_0
  if(length(x$a_0) <= 4){
  	print(tmp)
  }else{
  	tmp$cov <- diag(tmp$cov)
  }
}


#' Print the output of ci.drtmle
#' @export
#' @param x An object of class ci.drtmle
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print ci.drtmle	

print.ci.drtmle <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}

#' Print the output of ci.islptw
#' @export
#' @param x An object of class ci.islptw
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print ci.islptw 

print.ci.islptw <- function(x,digits = 3,...){
	tmp <- lapply(x, round, digits = digits)
	print(tmp)
}

#' Print the output of wald_test.drtmle
#' @export
#' @param x An object of class wald_test.drtmle
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print wald_test.drtmle

print.wald_test.drtmle <- function(x,digits = 3,...){
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
}

#' Print the output of wald_test.islptw
#' @export
#' @param x An object of class wald_test.islptw
#' @param digits Number of digits to round to 
#' @param ... Other options (not currently used)
#' @method print wald_test.islptw 

print.wald_test.islptw <- function(x,digits = 3,...){
  tmp <- lapply(x, round, digits = digits)
  print(tmp)
}


#--------------------------------------------------
# Plot methods
#--------------------------------------------------
#' Plot reduced dimension regression fits
#' 
#' @param x An object of class \code{"drtmle"}
#' @param nPoints Number of points to plot lines (increase for less bumpy plot,
#' decrease for faster evaluation).
#' @param a_0 For what value of a_0 should the plot be made for?
#' @param ... More arguments (not currently used) TO DO: Pass to plot? to lines?
#' @export
#' @method plot drtmle
#' @importFrom graphics axis lines par plot

plot.drtmle <- function(x, nPoints = 500, 
                        a_0 = x$a_0[1], ...){
  # ask to see next plot
  par(ask = TRUE)
  # check if returnModels is null
  if(is.null(x$QrnMod) & is.null(x$grnMod)){
    stop("Plot function only works if returnModels = TRUE.")
  }

  # which entry in x fits corresponds to this a_0
  listInd <- which(x$a_0 == a_0)

  # get local nuisance fits
  gn <- x$nuisance_drtmle$gnStar[[listInd]]
  Qn <- x$nuisance_drtmle$QnStar[[listInd]]

  #------------------
  # plot Qrn fit
  #------------------
  # number of fits (if no CV = 1, if CV > 1)
  nFit <- length(x$QrnMod)
  # xlim = range of gn 
  xl <- range(gn)
  # prediction points
  predP <- seq(xl[1],xl[2], length = nPoints)
  # get predictions back for each Qrn fit
  fit_Qrn <- lapply(x$QrnMod, function(y){
    newDat <- data.frame(gn = predP)
    if("glm" %in% class(y[[listInd]])){
      predict(y[[listInd]], newdata = newDat, type = "response")
    }else if(class(y[[listInd]]) == "SuperLearner"){
      pred <- predict(y[[listInd]], newdata = newDat)
      # get sl prediction if meta learning did not fail
      if(!all(y[[listInd]]$coef == 0)){
        pred$pred
      # otherwise get discrete super learner 
      }else{
        pred$library.predict[,which.min(y$cvRisk)]
      }
    }else{
      predict(y[[listInd]]$fit, newdata = newDat, type = "response")
    }
  })
  # get ylimits
  yl <- range(unlist(fit_Qrn))
  # set up empty plot
  plot(0, type = "n", xlim = xl, ylim = yl,
       xaxt = "n", yaxt = "n", bty = "n", 
       xlab = expression(g[n](W)), 
       ylab = expression("E[Y-"*Q[n](W)*" | "*g[n](W)*"]"))
  # add axes
  axis(side = 1); axis(side = 2)
  # add lines
  invisible(lapply(fit_Qrn, lines, x = predP, lwd = 2, col = "gray50"))

  #------------------
  # plot grn fit
  #------------------
  # only plot if univariate reduction
  reduction <- as.list(x$call)$reduction
  if(is.null(reduction)) reduction <- "univariate"
  if(reduction == "univariate"){
    # xlim = range of gn 
    xl <- range(Qn)
    # prediction points
    predP <- seq(xl[1], xl[2], length = nPoints)
    ## get fitted values of g_{n,r,1}
    fit_grn1 <- lapply(x$grnMod, function(y){
      newDat <- data.frame(Qn = predP)
      if("glm" %in% class(y[[listInd]]$fm1)){
        predict(y[[listInd]]$fm1, newdata = newDat, type = "response")
      }else if(class(y[[listInd]]$fm1) == "SuperLearner"){
        pred <- predict(y[[listInd]]$fm1, newdata = newDat)
        # get sl prediction if meta learning did not fail
        if(!all(y[[listInd]]$fm1$coef == 0)){
          pred$pred
        # otherwise get discrete super learner 
        }else{
          pred$library.predict[,which.min(y$cvRisk)]
        }
      }else{
        predict(y[[listInd]]$fm1$fit, newdata = newDat, type = "response")
      }
    })
    # get ylimits
    yl <- range(unlist(fit_grn1))
    # set up empty plot
    plot(0, type = "n", xlim = xl, ylim = yl,
         xaxt = "n", yaxt = "n", bty = "n", 
         xlab = expression(Q[n](W)), 
         ylab = expression("E[{"*A-g[n](W)*"} / "*g[n](W)*"} | "*Q[n](W)*"]"))
    # add axes
    axis(side = 1); axis(side = 2)
    # add lines
    invisible(lapply(fit_grn1, lines, x = predP, lwd = 2, col = "gray50"))

    ## get fitted values of g_{n,r,2}
    fit_grn2 <- lapply(x$grnMod, function(y){
      newDat <- data.frame(Qn = predP)
      if("glm" %in% class(y[[listInd]]$fm2)){
        predict(y[[listInd]]$fm2, newdata = newDat, type = "response")
      }else if(class(y[[listInd]]$fm2) == "SuperLearner"){
        pred <- predict(y[[listInd]]$fm2, newdata = newDat)
        # get sl prediction if meta learning did not fail
        if(!all(y[[listInd]]$fm2$coef == 0)){
          pred$pred
        # otherwise get discrete super learner 
        }else{
          pred$library.predict[,which.min(y$cvRisk)]
        }
      }else{
        predict(y[[listInd]]$fm2$fit, newdata = newDat, type = "response")
      }
    })
    # get ylimits
    yl <- range(unlist(fit_grn2))
    # set up empty plot
    plot(0, type = "n", xlim = xl, ylim = yl,
         xaxt = "n", yaxt = "n", bty = "n", 
         xlab = expression(Q[n](W)), 
         ylab = expression("E[A | "*Q[n](W)*"]"))
    # add axes
    axis(side = 1); axis(side = 2)
    # add lines
    invisible(lapply(fit_grn2, lines, x = predP, lwd = 2, col = "gray50"))
  }
}