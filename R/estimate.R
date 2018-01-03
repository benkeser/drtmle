#' estimateG
#'
#' Function to estimate propensity score
#'
#' @importFrom plyr alply
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or
#' 1)
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if
#' missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if
#' missing 1 if observed)
#' @param W A \code{data.frame} of named covariates
#' @param stratify A \code{boolean} indicating whether to estimate the missing
#' outcome regression separately
#' for observations with \code{A} equal to 0/1 (if \code{TRUE}) or to pool
#' across \code{A} (if \code{FALSE}).
#' @param SL_g A vector of characters describing the super learner library to be
#' used for each of the regression (\code{DeltaA}, \code{A}, and \code{DeltaY}).
#' To use the same regression for each of the regressions (or if there is no
#' missing data in \code{A} nor \code{Y}), a single library may be input.
#' @param tolg A numeric indicating the minimum value for estimates of the
#' propensity score.
#' @param verbose A boolean indicating whether to print status updates.
#' @param returnModels A boolean indicating whether to return model fits for the
#' outcome regression, propensity score, and reduced-dimension regressions.
#' @param glm_g A character describing a formula to be used in the call to
#' \code{glm} for the propensity score.
#' @param a_0 A vector of fixed treatment values at which to return marginal
#' mean estimates.
#' @param validRows A \code{list} of length \code{cvFolds} containing the row
#' indexes of observations to include in validation fold.
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#'

estimateG <- function(A, W, DeltaY, DeltaA, SL_g, glm_g, a_0, tolg,
                      stratify = FALSE, validRows = NULL, verbose = FALSE, returnModels = FALSE) {
  if (is.null(SL_g) & is.null(glm_g)) {
    stop("Specify Super Learner library or GLM formula for g")
  }
  if (!is.null(SL_g) & !is.null(glm_g)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      "Proceeding with Super Learner only."
    ))
    glm_g <- NULL
  }

  # subset data into training and validation sets
  if (length(validRows) != length(A)) {
    trainDeltaA <- DeltaA[-validRows]
    trainDeltaY <- DeltaY[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows, , drop = FALSE]
    validW <- W[validRows, , drop = FALSE]
    validA <- A[validRows]
    validDeltaA <- DeltaA[validRows]
    validDeltaY <- DeltaY[validRows]
  } else {
    trainA <- validA <- A
    trainW <- validW <- W
    trainDeltaA <- validDeltaA <- DeltaA
    trainDeltaY <- validDeltaY <- DeltaY
  }

  if (!is.null(SL_g)) {
    # check for names in SL_g
    namedSL_g <- c("DeltaA", "A", "DeltaY") %in% names(SL_g)
    # if none of the above names appear, then it is assumed that you want
    # to use SL_g for each of DeltaA, A, and Y
    if (!any(namedSL_g)) {
      SL_g <- list(DeltaA = SL_g, A = SL_g, DeltaY = SL_g)
    }
  } else if (!is.null(glm_g)) {
    namedglm_g <- c("DeltaA", "A", "DeltaY") %in% names(glm_g)
    # if none of the above names appear, then it is assumed that you want
    # to use glm_g for each of DeltaA, A, and Y
    if (!any(namedglm_g)) {
      glm_g <- list(DeltaA = glm_g, A = glm_g, DeltaY = glm_g)
    }
  }

  # -------------------------------
  # regression of DeltaA ~ W
  # -------------------------------
  # only fit this regression if there are some missing treatment assignments
  if (!all(DeltaA == 1)) {
    # if super learner library is specified, fit a super learner
    if (!is.null(SL_g)) {
      # if the SL_g$DeltaA is of length > 1, then call SuperLearner
      if (length(SL_g$DeltaA) > 1 | is.list(SL_g$DeltaA)) {
        fm_DeltaA <- SuperLearner::SuperLearner(
          Y = trainDeltaA,
          X = trainW, newX = validW, family = stats::binomial(),
          SL.library = SL_g$DeltaA, verbose = verbose,
          method = "method.CC_nloglik_mod"
        )
        # get predicted probability of missing treatment
        gn_DeltaA <- fm_DeltaA$SL.predict
      } else if (!is.list(SL_g$DeltaA) & length(SL_g$DeltaA) == 1) {
        fm_DeltaA <- do.call(SL_g$DeltaA, args = list(
          Y = trainDeltaA,
          X = trainW, newX = validW,
          obsWeights = rep(1, length(trainA)),
          family = stats::binomial()
        ))
        gn_DeltaA <- fm_DeltaA$pred
      }
    } # end if SuperLearner loop
    if (!is.null(glm_g)) {
      thisDat <- data.frame(DeltaA = trainDeltaA, trainW)
      fm_DeltaA <- stats::glm(stats::as.formula(paste0(
        "DeltaA~",
        glm_g$DeltaA
      )), data = thisDat, family = stats::binomial())
      gn_DeltaA <- stats::predict(
        fm_DeltaA, type = "response",
        newdata = data.frame(DeltaA = validDeltaA, validW)
      )
    }
    # name for returned models
    name_DeltaA <- "DeltaA ~ W"
  } else {
    # if all DeltaA==1 then put NULL model and 1 predictions
    fm_DeltaA <- NULL
    name_DeltaA <- ""
    gn_DeltaA <- rep(1, length(validDeltaA))
  }

  # -----------------------------------
  # fitting A ~ W | DeltaA = 1
  # -----------------------------------
  # if a super learner library is specified, fit the super learner
  if (!is.null(SL_g)) {
    # if the library is of length > 1, then call SuperLearner
    if (length(SL_g$A) > 1 | is.list(SL_g$A)) {
      # if there are only two unique values of A, then only need one fit
      if (length(a_0) == length(unique(A)) &
        length(unique(A[!is.na(A)])) == 2) {
        fm_A <- SuperLearner::SuperLearner(
          Y = as.numeric(trainA[trainDeltaA == 1] == a_0[1]),
          X = trainW[trainDeltaA == 1, , drop = FALSE], newX = validW,
          family = stats::binomial(), SL.library = SL_g$A,
          verbose = verbose, method = "method.CC_nloglik_mod"
        )
        gn_A <- vector(mode = "list", length = 2)
        gn_A[[1]] <- fm_A$SL.predict
        gn_A[[2]] <- 1 - gn_A[[1]]
        # name for this model
        name_A <- paste0("I(A = ", a_0[1], ") ~ W | DeltaA == 1")
        # if there are more than two unique values of A, then we need
        # more than one call to super learner
      } else {
        a_ct <- 0
        gn_A <- vector(mode = "list", length = length(a_0))
        fm_A <- vector(mode = "list", length = length(a_0) - 1)
        name_A <- rep(NA, length(a_0) - 1)
        for (a in a_0[1:(length(a_0) - 1)]) {
          # determine who to include in the regression for this outcome
          if (a_ct == 0) {
            include <- rep(TRUE, length(trainA))
          } else {
            include <- !(trainA %in% a_0[1:a_ct])
          }
          # now exclude people with DeltaA = 0
          include[trainDeltaA == 0] <- FALSE
          # fit super learner
          tmp_fm <- SuperLearner::SuperLearner(
            Y = as.numeric(trainA[include] == a),
            X = trainW[include, , drop = FALSE], newX = validW,
            family = stats::binomial(), SL.library = SL_g$A,
            verbose = verbose
          )
          # get predictions
          tmp_pred <- tmp_fm$SL.pred
          if (a_ct != 0) {
            gn_A[[a_ct + 1]] <- tmp_pred * Reduce(
              "*",
              lapply(gn_A[1:a_ct], function(x) {
                1 - x
              })
            )
          } else {
            gn_A[[a_ct + 1]] <- tmp_pred
          }
          fm_A[[a_ct + 1]] <- tmp_fm
          name_A[a_ct + 1] <- paste0("I(A = ", a, ") ~ W | DeltaA == 1")
          a_ct <- a_ct + 1
        }
        # add in final predictions
        gn_A[[a_ct + 1]] <- 1 - Reduce("+", gn_A[1:a_ct])
      }
    } else if (!is.list(SL_g$A) & length(SL_g$A) == 1) {
      if (length(a_0) == length(unique(A[!is.na(A)])) &
        length(unique(A[!is.na(A)])) == 2) {
        gn_A <- vector(mode = "list", length = 2)
        fm_A <- do.call(SL_g$A, args = list(
          Y = as.numeric(
            trainA[trainDeltaA == 1] == a_0[1]
          ),
          X = trainW[trainDeltaA == 1, , drop = FALSE], newX = validW,
          obsWeights = rep(1, length(trainA[trainDeltaA == 1])),
          family = stats::binomial()
        ))
        gn_A[[1]] <- fm_A$pred
        gn_A[[2]] <- 1 - fm_A$pred
        name_A <- paste0("I(A = ", a_0[1], ") ~ W | DeltaA == 1")
      } else {
        a_ct <- 0
        gn_A <- vector(mode = "list", length = length(a_0))
        fm_A <- vector(mode = "list", length = length(a_0) - 1)
        name_A <- rep(NA, length(a_0) - 1)
        for (a in a_0[1:(length(a_0) - 1)]) {
          # determine who to include in the regression for this outcome
          if (a_ct == 0) {
            include <- rep(TRUE, length(trainA))
          } else {
            include <- !(trainA %in% a_0[1:a_ct])
          }
          # set missing treatment people to FALSE
          include[trainDeltaA == 0] <- FALSE
          # fit super learner
          tmp_fm <- do.call(SL_g$A, args = list(
            Y = as.numeric(
              trainA[include] == a
            ), X = trainW[include, , drop = FALSE],
            newX = validW, obsWeights = rep(1, length(trainA[include])),
            family = stats::binomial()
          ))
          # get predictions
          tmp_pred <- tmp_fm$pred
          if (a_ct != 0) {
            gn_A[[a_ct + 1]] <- tmp_pred * Reduce(
              "*",
              lapply(gn_A[1:a_ct], function(x) {
                1 - x
              })
            )
          } else {
            gn_A[[a_ct + 1]] <- tmp_pred
          }
          fm_A[[a_ct + 1]] <- tmp_fm
          name_A[a_ct + 1] <- paste0("I(A = ", a, ") ~ W | DeltaA == 1")
          a_ct <- a_ct + 1
        }
        # add in final predictions
        gn_A[[a_ct + 1]] <- 1 - Reduce("+", gn_A[1:a_ct])
      }
    }
  }
  # ----------------------------------------------------------------------
  # GLM
  # ----------------------------------------------------------------------
  if (!is.null(glm_g)) {
    if (length(a_0) == length(unique(A)) &
      length(unique(A[!is.na(A)])) == 2) {
      thisDat <- data.frame(A = as.numeric(trainA[trainDeltaA == 1] ==
        a_0[1]), trainW[trainDeltaA == 1, , drop = FALSE])
      fm_A <- stats::glm(
        stats::as.formula(paste0("A~", glm_g$A)),
        data = thisDat, family = stats::binomial()
      )
      gn_A <- vector(mode = "list", length = 2)
      name_A <- paste0("I(A = ", a_0[1], ") ~ W | DeltaA == 1")
      gn_A[[1]] <- stats::predict(fm_A, newdata = data.frame(
        A = validA, validW
      ), type = "response")
      gn_A[[2]] <- 1 - gn_A[[1]]
    } else {
      a_ct <- 0
      gn_A <- vector(mode = "list", length = length(a_0))
      fm_A <- vector(mode = "list", length = length(a_0) - 1)
      name_A <- rep(NA, length(a_0) - 1)
      for (a in a_0[1:(length(a_0) - 1)]) {
        # determine who to include in the regression for this outcome
        if (a_ct == 0) {
          include <- rep(TRUE, length(A))
        } else {
          include <- !(A %in% a_0[1:a_ct])
        }
        # don't include folks with missing treatment
        include[trainDeltaA == 0] <- FALSE
        # fit super learner
        thisDat <- data.frame(
          as.numeric(trainA[include] == a),
          trainW[include, , drop = FALSE]
        )
        colnames(thisDat) <- c("A", colnames(W))
        tmp_fm <- stats::glm(
          stats::as.formula(paste0("A~", glm_g)),
          data = thisDat, family = stats::binomial()
        )
        tmp_pred <- stats::predict(tmp_fm, newdata = data.frame(
          A = validA, validW
        ), type = "response")
        # get predictions
        if (a_ct != 0) {
          gn_A[[a_ct + 1]] <- tmp_pred * Reduce(
            "*",
            lapply(gn_A[1:a_ct], function(x) {
              1 - x
            })
          )
        } else {
          gn_A[[a_ct + 1]] <- tmp_pred
        }
        fm_A[[a_ct + 1]] <- tmp_fm
        name_A[a_ct + 1] <- paste0("I(A = ", a, ") ~ W | DeltaA == 1")
        a_ct <- a_ct + 1
      } # end for loop over treatment levels
      # add in final predictions
      gn_A[[a_ct + 1]] <- 1 - Reduce("+", gn_A[1:a_ct])
    } # end multi-level treatment if
  } # end glm_g if

  # -------------------------------------
  # fit DeltaY ~ W + A | DeltaA = 1
  # -------------------------------------
  # only fit this regression if there are some missing outcomes
  if (!all(DeltaY == 1)) {
    # only include people with DeltaA == 1
    include <- (trainDeltaA == 1)

    # if super learner library is specified, fit a super learner
    if (!is.null(SL_g)) {
      # if the SL_g$DeltaY is of length > 1, then call SuperLearner
      if (length(SL_g$DeltaY) > 1 | is.list(SL_g$DeltaY)) {
        # if no stratify, then fit DeltaY ~ W | DeltaA = 1 in each
        # level of A
        if (stratify) {
          fm_DeltaY <- vector(mode = "list", length = length(a_0))
          gn_DeltaY <- vector(mode = "list", length = length(a_0))
          name_DeltaY <- rep(NA, length(a_0))
          a_ct <- 0
          for (a in a_0) {
            a_ct <- a_ct + 1

            # only include people with A == a and DeltaA == 1
            include2 <- (trainA == a)
            include2[is.na(include2)] <- FALSE

            # fit super learner
            fm_DeltaY[[a_ct]] <- SuperLearner::SuperLearner(
              Y = trainDeltaY[include & include2],
              X = trainW[include & include2, , drop = FALSE],
              newX = validW, family = stats::binomial(),
              SL.library = SL_g$DeltaY, verbose = verbose,
              method = "method.CC_nloglik_mod"
            )
            # name the fit
            name_DeltaY[a_ct] <- paste0(
              "DeltaY ~ W | DeltaA == 1",
              " & A == ", a
            )
            # get predictions back on everybody
            gn_DeltaY[[a_ct]] <- fm_DeltaY[[a_ct]]$SL.predict
          } # end loop over treatment levels
          # if not stratified, fit a single regression pooling over
          # levels of A
        } else {
          # fit super learner
          fm_DeltaY <- SuperLearner::SuperLearner(
            Y = trainDeltaY[include], X = data.frame(
              A = trainA[include], trainW[include, , drop = FALSE]
            ),
            family = stats::binomial(), SL.library = SL_g$DeltaY,
            verbose = verbose, method = "method.CC_nloglik_mod"
          )

          # get predictions back setting A = a for every a in a_0
          gn_DeltaY <- vector(mode = "list", length = length(a_0))
          name_DeltaY <- paste0("DeltaY ~ W + A | DeltaA == 1")
          a_ct <- 0
          for (a in a_0) {
            a_ct <- a_ct + 1
            gn_DeltaY[[a_ct]] <- stats::predict(
              fm_DeltaY,
              onlySL = TRUE, newdata = data.frame(A = a, validW)
            )$pred
          }
        } # end if !stratify
        # if SL_g$DeltaY only a single algorithm, then call directly
      } else if (!is.list(SL_g$DeltaY) & length(SL_g$DeltaY) == 1) {
        # if no stratify, then fit DeltaY ~ W | DeltaA = 1 in
        # each level of A
        if (stratify) {
          fm_DeltaY <- vector(mode = "list", length = length(a_0))
          gn_DeltaY <- vector(mode = "list", length = length(a_0))
          name_DeltaY <- rep(NA, length(a_0))
          a_ct <- 0
          for (a in a_0) {
            a_ct <- a_ct + 1

            # only include people with A == a
            include2 <- (trainA == a)
            include2[is.na(include2)] <- FALSE
            # make call to algorithm
            fm_DeltaY[[a_ct]] <- do.call(SL_g$DeltaY, args = list(
              Y = trainDeltaY[include & include2],
              X = trainW[include & include2, , drop = FALSE],
              newX = validW, obsWeights = rep(
                1,
                length(trainA[include & include2])
              ),
              family = stats::binomial()
            ))
            name_DeltaY[a_ct] <- paste0(
              "DeltaY ~ W | DeltaA == 1",
              " & A == ", a
            )
            # get predictions
            gn_DeltaY[[a_ct]] <- fm_DeltaY[[a_ct]]$pred
          } # end loop over a_0
        } else {
          # end if stratify call algorithm to fit pooled estimate
          fm_DeltaY <- do.call(SL_g$DeltaY, args = list(
            Y = trainDeltaY[include],
            X = data.frame(
              A = trainA[include],
              trainW[include, , drop = FALSE]
            ), newX = data.frame(
              A = validA, validW
            ), obsWeights = rep(
              1,
              length(trainA[include])
            ), family = stats::binomial()
          ))
          name_DeltaY <- paste0("DeltaY ~ W + A | DeltaA == 1")
          # loop to get predictions setting A = a
          gn_DeltaY <- vector(mode = "list", length = length(a_0))
          a_ct <- 0
          for (a in a_0) {
            a_ct <- a_ct + 1
            gn_DeltaY[[a_ct]] <- stats::predict(
              fm_DeltaY$fit,
              newdata = data.frame(A = a, validW)
            )
          }
        } # end !stratify
      } # end if one algorithm loop
    } # end SL_g not null if
    if (!is.null(glm_g)) {
      if (stratify) {
        fm_DeltaY <- vector(mode = "list", length = length(a_0))
        gn_DeltaY <- vector(mode = "list", length = length(a_0))
        name_DeltaY <- rep(NA, length = length(a_0))
        a_ct <- 0
        for (a in a_0) {
          a_ct <- a_ct + 1

          # only include people with A == a and DeltaA == 1
          include2 <- (trainA == a)
          include2[is.na(include2)] <- FALSE
          fm_DeltaY[[a_ct]] <- stats::glm(stats::as.formula(
            paste0(
              "trainDeltaY[include & include2]~",
              glm_g$DeltaY
            )
          ), data = data.frame(trainW[include &
            include2, , drop = FALSE]), family = stats::binomial())
          name_DeltaY[a_ct] <- paste0(
            "DeltaY ~ W | DeltaA == 1 ",
            "& A == ", a
          )
          # get predictions back for everyone
          gn_DeltaY[[a_ct]] <- stats::predict(
            fm_DeltaY[[a_ct]],
            newdata = validW, type = "response"
          )
        } # end loop over treatments
      } else {
        # end stratified glm fit glm in everyone with DeltaA == 1
        fm_DeltaY <- stats::glm(
          stats::as.formula(paste0(
            "trainDeltaY[include]~", glm_g$DeltaY
          )),
          data = data.frame(A = trainA[include], trainW[
            include,
            , drop = FALSE
          ]), family = stats::binomial()
        )
        name_DeltaY <- paste0("DeltaY ~ W + A | DeltaA == 1")
        # get predictions back setting A = a
        gn_DeltaY <- vector(mode = "list", length = length(a_0))
        a_ct <- 0
        for (a in a_0) {
          a_ct <- a_ct + 1
          gn_DeltaY[[a_ct]] <- stats::predict(
            fm_DeltaY,
            newdata = data.frame(A = a, validW), type = "response"
          )
        } # end loop over treatments
      } # end !stratified glm
    } # end glm if
  } else {
    # if all DeltaY==1 then NULL model and 1 pred.
    fm_DeltaY <- NULL
    name_DeltaY <- ""
    gn_DeltaY <- vector(mode = "list", length = length(a_0))
    for (i in 1:length(a_0)) {
      gn_DeltaY[[i]] <- rep(1, length(validDeltaY))
    }
  }

  # ------------------------------------------------------
  # combine estimates into a single propensity score
  # ------------------------------------------------------
  gn <- mapply(gn_A = gn_A, gn_DeltaY = gn_DeltaY, FUN = function(gn_A,
                                                                  gn_DeltaY) {
    gn_A * gn_DeltaY * gn_DeltaA
  }, SIMPLIFY = FALSE)

  # truncate too-small predictions
  gn <- lapply(gn, function(g) {
    g[g < tolg] <- tolg
    g
  })

  out <- list(est = gn, fm = NULL)
  if (returnModels) {
    names(fm_A) <- name_A
    if (!is.null(fm_DeltaA)) {
      names(fm_DeltaA) <- name_DeltaA
    }
    if (!is.null(fm_DeltaY)) {
      names(fm_DeltaY) <- name_DeltaY
    }
    out$fm <- list(DeltaA = fm_DeltaA, A = fm_A, DeltaY = fm_DeltaY)
  }
  return(out)
}



#' estimateQ
#'
#' Function to estimate initial outcome regression
#'
#' @param Y A vector of continuous or binary outcomes.
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or
#' 1).
#' @param W A \code{data.frame} of named covariates.
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if
#' missing 1 if observed).
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if
#' missing 1 if observed).
#' @param SL_Q A vector of characters or a list describing the Super Learner
#' library to be used for the outcome regression.
#' @param verbose A boolean indicating whether to print status updates.
#' @param returnModels A boolean indicating whether to return model fits for the
#' outcome regression, propensity score, and reduced-dimension regressions.
#' @param glm_Q A character describing a formula to be used in the call to
#' \code{glm} for the outcome regression.
#' @param a_0 A list of fixed treatment values
#' @param family A character passed to \code{SuperLearner}
#' @param stratify A \code{boolean} indicating whether to estimate the outcome
#' regression separately for observations with \code{A} equal to 0/1 (if
#' \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param validRows A \code{list} of length \code{cvFolds} containing the row
#' indexes of observations to include in validation fold.
#' @param ... Additional arguments (not currently used)
#'
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#'
estimateQ <- function(Y, A, W, DeltaA, DeltaY, SL_Q, glm_Q, a_0, stratify,
                      family, verbose = FALSE, returnModels = FALSE, validRows = NULL,
                      ...) {
  if (is.null(SL_Q) & is.null(glm_Q)) {
    stop("Specify Super Learner library or GLM formula for Q")
  }
  if (!is.null(SL_Q) & !is.null(glm_Q)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      " Proceeding with Super Learner only."
    ))
    glm_Q <- NULL
  }
  # subset data into training and validation sets
  if (length(validRows) != length(Y)) {
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows, , drop = FALSE]
    trainDeltaA <- DeltaA[-validRows]
    trainDeltaY <- DeltaY[-validRows]
    validW <- W[validRows, , drop = FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
    validDeltaY <- DeltaY[validRows]
    validDeltaA <- DeltaA[validRows]
  } else {
    trainA <- validA <- A
    trainW <- validW <- W
    trainY <- validY <- Y
    trainDeltaA <- validDeltaA <- DeltaA
    trainDeltaY <- validDeltaY <- DeltaY
  }

  # include only DeltaA = 1 and DeltaY = 1 folks
  include <- (trainDeltaA == 1) & (trainDeltaY == 1)

  # Super Learner
  if (!is.null(SL_Q)) {
    if (!stratify) {
      if (length(SL_Q) > 1 | is.list(SL_Q)) {
        fm <- SuperLearner::SuperLearner(
          Y = trainY[include],
          X = data.frame(A = trainA, trainW)[include, , drop = FALSE],
          verbose = verbose, family = family, SL.library = SL_Q,
          method = ifelse(family$family == "binomial",
            "method.CC_nloglik_mod", "method.CC_LS_mod"
          )
        )

        Qn <- alply(a_0, 1, function(x) {
          stats::predict(
            fm, newdata = data.frame(A = x, validW),
            onlySL = TRUE
          )[[1]]
        })
      } else if (length(SL_Q) == 1) {
        fm <- do.call(SL_Q, args = list(
          Y = trainY[include],
          X = data.frame(A = trainA, trainW)[include, , drop = FALSE],
          verbose = verbose, newX = data.frame(A = validA, validW),
          obsWeights = rep(1, length(trainA[include])),
          family = family
        ))
        Qn <- alply(a_0, 1, function(x) {
          stats::predict(object = fm$fit, newdata = data.frame(
            A = x,
            validW
          ))
        })
      }
    } else {
      if (length(SL_Q) > 1 | is.list(SL_Q)) {
        tmp <- plyr::alply(a_0, 1, function(x) {
          include2 <- trainA == x
          # handle NAs properly
          include2[is.na(include2)] <- FALSE
          fm <- SuperLearner::SuperLearner(
            Y = trainY[include2 & include],
            X = trainW[include2 & include, , drop = FALSE],
            newX = validW, verbose = verbose, family = family,
            SL.library = SL_Q, method = ifelse(family$family ==
              "binomial", "method.CC_nloglik_mod", "method.CC_LS_mod")
          )
          list(est = fm$SL.predict, fm = fm)
        })
        Qn <- lapply(tmp, "[[", 1)
        fm <- lapply(tmp, "[[", 2)
      } else if (length(SL_Q) == 1) {
        tmp <- plyr::alply(a_0, 1, function(x) {
          include2 <- trainA == x
          # handle NAs properly
          include2[is.na(include2)] <- FALSE
          # call function
          fm <- do.call(SL_Q, args = list(
            Y = trainY[include2 & include],
            X = trainW[include2 & include, , drop = FALSE],
            newX = validW, verbose = verbose,
            obsWeights = rep(1, sum(include2 & include)),
            family = family
          ))
          list(est = fm$pred, fm = fm)
        })
        Qn <- lapply(tmp, "[[", 1)
        fm <- lapply(tmp, "[[", 2)
      }
    }
  }

  # GLM
  if (!is.null(glm_Q)) {
    if (!stratify) {
      fm <- stats::glm(
        stats::as.formula(paste0("Y~", glm_Q)),
        data = data.frame(Y = trainY, A = trainA, trainW)[
          include, ,
          drop = FALSE
        ], family = family
      )
      Qn <- plyr::alply(matrix(a_0), 1, function(a, fm) {
        stats::predict(
          fm, newdata = data.frame(A = a, validW),
          type = "response"
        )
      }, fm = fm)
    } else {
      tmp <- plyr::alply(matrix(a_0), 1, function(a) {
        include2 <- trainA == a
        # handle NAs properly
        include2[is.na(include2)] <- FALSE
        fm <- stats::glm(
          stats::as.formula(paste0(
            "trainY[include2 & include] ~ ", glm_Q
          )),
          data = trainW[include2 & include, , drop = FALSE],
          family = family
        )
        return(list(est = stats::predict(
          fm, newdata = validW,
          type = "response"
        ), fm = fm))
      })
      Qn <- lapply(tmp, "[[", 1)
      fm <- lapply(tmp, "[[", 2)
    }
  }
  out <- list(est = Qn, fm = NULL)
  if (returnModels) {
    out$fm <- fm
  }
  return(out)
}

#' estimateQrn
#'
#' Estimates the reduced dimension regressions necessary for the
#' fluctuations of g
#'
#'
#' @param Y A vector of continuous or binary outcomes.
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or
#'  1)
#' @param W A \code{data.frame} of named covariates
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if
#' missing 1 if observed)
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if
#' missing 1 if observed)
#' @param Qn A list of outcome regression estimates evaluated on observed data.
#' If NULL then 0 is used for all Qn (as is needed to estimate reduced dimension
#' regression for adaptive_iptw)
#' @param gn A list of propensity regression estimates evaluated on observed
#' data
#' @param SL_Qr A vector of characters or a list describing the Super Learner
#' library to be used for the first reduced-dimension regression.
#' @param glm_Qr A character describing a formula to be used in the call to
#' \code{glm} for the first reduced-dimension regression. Ignored if
#' \code{SL_gr!=NULL}.
#' @param a_0 A list of fixed treatment values.
#' @param returnModels A boolean indicating whether to return model fits for the
#' outcome regression, propensity score, and reduced-dimension regressions.
#' @param family Should be gaussian() unless called from adaptive_iptw with
#' binary \code{Y}.
#' @param validRows A \code{list} of length \code{cvFolds} containing the row
#' indexes of observations to include in validation fold.
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula gaussian binomial

estimateQrn <- function(Y, A, W, DeltaA, DeltaY, Qn, gn, glm_Qr, SL_Qr,
                        family = stats::gaussian(), a_0, returnModels, validRows = NULL) {

  # if estimateQrn is called in adaptive_iptw, then Qn will enter as NULL.
  # Here we fill its value to 0 so that we estimate the correct nuisance
  # parameter for adaptive_iptw
  if (is.null(Qn)) {
    Qn <- vector(mode = "list", length = length(a_0))
    for (i in 1:length(a_0)) {
      Qn[[i]] <- rep(0, length(Y))
    }
  }

  # subset data into training and validation sets
  if (length(validRows) != length(Y)) {
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows, , drop = FALSE]
    trainDeltaA <- DeltaA[-validRows]
    trainDeltaY <- DeltaY[-validRows]
    train_gn <- lapply(gn, "[", -validRows)
    train_Qn <- lapply(Qn, "[", -validRows)
    validW <- W[validRows, , drop = FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
    validDeltaA <- DeltaA[-validRows]
    validDeltaY <- DeltaY[-validRows]
    valid_gn <- lapply(gn, "[", validRows)
    valid_Qn <- lapply(Qn, "[", validRows)
  } else {
    trainA <- validA <- A
    trainW <- validW <- W
    trainY <- validY <- Y
    trainDeltaY <- validDeltaY <- DeltaY
    trainDeltaA <- validDeltaA <- DeltaA
    train_gn <- valid_gn <- gn
    train_Qn <- valid_Qn <- Qn
  }

  if (is.null(SL_Qr) & is.null(glm_Qr)) {
    stop("Specify Super Learner library or GLM formula for Qr")
  }
  if (!is.null(SL_Qr) & !is.null(glm_Qr)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      "Proceeding with Super Learner only."
    ))
    glm_Qr <- NULL
  }
  # Super Learner
  if (!is.null(SL_Qr)) {
    Qrn <- mapply(
      a = a_0, train_g = train_gn, train_Q = train_Qn,
      valid_g = valid_gn, valid_Q = valid_Qn, SIMPLIFY = FALSE,
      FUN = function(a, train_g, train_Q, valid_g, valid_Q) {
        Aeqa <- trainA == a
        Aeqa[is.na(Aeqa)] <- FALSE
        if (length(unique(train_g)) == 1) {
          warning(paste0(
            "Only one unique value of gn", a,
            ". Using empirical average as Qr estimate."
          ))
          m1 <- mean((trainY - train_Q)[Aeqa & trainDeltaA == 1 &
            trainDeltaY == 1])
          est <- rep(m1, length(validY))
          fm <- list(fit = list(object = m1), pred = NULL)
          class(fm$fit) <- "SL.mean"
        } else {
          if (length(SL_Qr) > 1) {
            suppressWarnings(fm <- SuperLearner::SuperLearner(
              Y = (trainY - train_Q)[Aeqa & trainDeltaA == 1 &
                trainDeltaY == 1], X = data.frame(gn = train_g[Aeqa &
                trainDeltaA == 1 & trainDeltaY == 1]), newX = data.frame(
                gn = valid_g
              ), family = family, SL.library = SL_Qr,
              method = "method.CC_LS_mod"
            ))
            # if all weights = 0, use discrete SL
            if (!all(fm$coef == 0)) {
              est <- fm$SL.predict
            } else {
              est <- fm$library.predict[, which.min(fm$cvRisk)]
            }
          } else if (length(SL_Qr) == 1) {
            fm <- do.call(SL_Qr, args = list(
              Y = (trainY - train_Q)[Aeqa & trainDeltaA == 1 &
                trainDeltaY == 1], X = data.frame(gn = train_g[Aeqa &
                trainDeltaA == 1 & trainDeltaY == 1]),
              family = family, newX = data.frame(gn = valid_g),
              obsWeights = rep(1, length(trainY[Aeqa &
                trainDeltaA == 1 & trainDeltaY == 1]))
            ))
            est <- fm$pred
          }
        }
        out <- list(est = est, fm = NULL)
        if (returnModels) {
          out$fm <- fm
        }
        return(out)
      }
    )
  }

  # GLM
  if (!is.null(glm_Qr)) {
    Qrn <- mapply(
      a = a_0, train_g = train_gn, train_Q = train_Qn,
      valid_g = valid_gn, valid_Q = valid_Qn, SIMPLIFY = FALSE,
      FUN = function(a, train_g, train_Q, valid_g, valid_Q) {
        Aeqa <- trainA == a
        Aeqa[is.na(Aeqa)] <- FALSE
        if (length(unique(train_g)) == 1) {
          warning(paste0(
            "Only one unique value of gn", a,
            ". Using empirical average as Qr estimate."
          ))
          glm_Qr <- "1"
        }
        fm <- stats::glm(
          stats::as.formula(paste0("Qrn ~", glm_Qr)),
          data = data.frame(
            Qrn = (trainY - train_Q)[Aeqa &
              trainDeltaY == 1 & trainDeltaY == 1],
            gn = train_g[Aeqa & trainDeltaY == 1 & trainDeltaY == 1]
          ),
          family = family
        )
        est <- stats::predict(
          fm, newdata = data.frame(gn = valid_g),
          type = "response"
        )
        out <- list(est = est, fm = NULL)
        if (returnModels) {
          out$fm <- fm
        }
        return(out)
      }
    )
  }
  # return estimates and models
  return(list(est = lapply(Qrn, function(x) {
    x$est
  }), fm = lapply(Qrn, function(x) {
    fm <- x$fm
  })))
  Qrn
}

#' estimategrn
#'
#' Estimates the reduced dimension regressions necessary for the additional
#' fluctuations.
#'
#' @param Y A vector of continuous or binary outcomes.
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or
#'  1).
#' @param W A \code{data.frame} of named covariates.
#' @param DeltaY Indicator of missing outcome (assumed to be equal to 0 if
#' missing 1 if observed).
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if
#' missing 1 if observed).
#' @param Qn A list of outcome regression estimates evaluated on observed data.
#' @param gn A list of propensity regression estimates evaluated on observed
#' data.
#' @param SL_gr A vector of characters or a list describing the Super Learner
#' library to be used for the reduced-dimension propensity score.
#' @param glm_gr A character describing a formula to be used in the call to
#' \code{glm} for the second reduced-dimension regression. Ignored if
#' \code{SL_gr!=NULL}.
#' @param reduction A character equal to \code{'univariate'} for a univariate
#' misspecification correction or \code{'bivariate'} for the bivariate version.
#' @param tolg A numeric indicating the minimum value for estimates of the
#' propensity score.
#' @param a_0 A list of fixed treatment values .
#' @param returnModels A boolean indicating whether to return model fits for the
#' outcome regression, propensity score, and reduced-dimension regressions.
#' @param validRows A \code{list} of length \code{cvFolds} containing the row
#' indexes of observations to include in validation fold.
#'
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula

estimategrn <- function(Y, A, W, DeltaA, DeltaY, Qn, gn, SL_gr, tolg, glm_gr,
                        a_0, reduction, returnModels, validRows) {
  if (length(validRows) != length(Y)) {
    trainY <- Y[-validRows]
    trainA <- A[-validRows]
    trainW <- W[-validRows, , drop = FALSE]
    trainDeltaA <- DeltaA[-validRows]
    trainDeltaY <- DeltaY[-validRows]
    train_gn <- lapply(gn, "[", -validRows)
    train_Qn <- lapply(Qn, "[", -validRows)
    validW <- W[validRows, , drop = FALSE]
    validA <- A[validRows]
    validY <- Y[validRows]
    validDeltaA <- DeltaA[-validRows]
    validDeltaY <- DeltaY[-validRows]
    valid_gn <- lapply(gn, "[", validRows)
    valid_Qn <- lapply(Qn, "[", validRows)
  } else {
    trainA <- validA <- A
    trainW <- validW <- W
    trainY <- validY <- Y
    trainDeltaY <- validDeltaY <- DeltaY
    trainDeltaA <- validDeltaA <- DeltaA
    train_gn <- valid_gn <- gn
    train_Qn <- valid_Qn <- Qn
  }

  if (is.null(SL_gr) & is.null(glm_gr)) {
    stop("Specify Super Learner library or GLM formula for gr")
  }
  if (!is.null(SL_gr) & !is.null(glm_gr)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      "Proceeding with Super Learner only."
    ))
    glm_gr <- NULL
  }
  # Super Learner
  if (!is.null(SL_gr)) {
    grn <- mapply(
      a = a_0, train_Q = train_Qn, train_g = train_gn,
      valid_Q = valid_Qn, valid_g = valid_gn, SIMPLIFY = FALSE,
      FUN = function(a, train_Q, train_g, valid_Q, valid_g) {
        Aeqa <- trainA == a
        Aeqa[is.na(Aeqa)] <- FALSE
        if (length(unique(train_Q)) == 1) {
          warning(paste0(
            "Only one unique value of Qn.",
            "Proceeding with empirical mean for grn"
          ))
          if (reduction == "univariate") {
            m1 <- mean((as.numeric(Aeqa & trainDeltaA == 1 &
              trainDeltaY == 1) - train_g) / train_g)
            grn1 <- rep(m1, length(validY))
            m2 <- mean(as.numeric(Aeqa & trainDeltaA == 1 &
              trainDeltaY == 1))
            grn2 <- rep(m2, length(validY))
            grn2[grn2 < tolg] <- tolg
            fm1 <- list(fit = list(object = m1), pred = NULL)
            class(fm1$fit) <- "SL.mean"
            fm2 <- list(fit = list(object = m2), pred = NULL)
            class(fm2$fit) <- "SL.mean"
          } else if (reduction == "bivariate") {
            m2 <- mean(as.numeric(Aeqa & trainDeltaA == 1 &
              trainDeltaY == 1))
            grn2 <- rep(m2, length(validY))
            grn2[grn2 < tolg] <- tolg
            fm2 <- list(fit = list(object = m2), pred = NULL)
            class(fm2$fit) <- "SL.mean"
            fm1 <- NULL
            grn1 <- rep(NA, length(validY))
          }
        } else {
          if (length(SL_gr) > 1) {
            if (reduction == "univariate") {
              fm1 <- SuperLearner::SuperLearner(
                Y = (as.numeric(Aeqa & trainDeltaA == 1 &
                  trainDeltaY == 1) - train_g) / train_g,
                X = data.frame(Qn = train_Q),
                newX = data.frame(Qn = valid_Q),
                family = stats::gaussian(), SL.library = SL_gr,
                method = "method.CC_LS_mod"
              )
              fm2 <- SuperLearner::SuperLearner(
                Y = as.numeric(Aeqa &
                  trainDeltaA == 1 & trainDeltaY == 1),
                X = data.frame(Qn = train_Q),
                newX = data.frame(Qn = valid_Q),
                family = stats::binomial(), SL.library = SL_gr,
                method = "method.CC_nloglik_mod"
              )
              if (!all(fm1$coef == 0)) {
                grn1 <- fm1$SL.predict
              } else {
                grn1 <- fm1$library.predict[, which.min(fm1$cvRisk)]
              }

              if (!all(fm2$coef == 0)) {
                grn2 <- fm2$SL.predict
              } else {
                grn2 <- fm2$library.predict[, which.min(fm2$cvRisk)]
              }
              grn2[grn2 < tolg] <- tolg
            } else if (reduction == "bivariate") {
              fm2 <- SuperLearner::SuperLearner(
                Y = as.numeric(Aeqa &
                  trainDeltaA == 1 & trainDeltaY == 1),
                X = data.frame(Qn = train_Q, gn = train_g),
                newX = data.frame(Qn = valid_Q, gn = valid_g),
                family = stats::binomial(),
                SL.library = SL_gr, method = "method.CC_nloglik_mod"
              )
              if (!all(fm2$coef == 0)) {
                grn2 <- fm2$SL.predict
              } else {
                grn2 <- fm2$library.predict[, which.min(fm2$cvRisk)]
              }
              grn2[grn2 < tolg] <- tolg
              fm1 <- NULL
              grn1 <- rep(NA, length(validY))
            }
          } else if (length(SL_gr) == 1) {
            if (reduction == "univariate") {
              fm1 <- do.call(SL_gr, args = list(
                Y = (as.numeric(Aeqa &
                  trainDeltaA == 1 & trainDeltaY == 1) - train_g) / train_g,
                X = data.frame(Qn = train_Q), obsWeights = rep(
                  1,
                  length(trainA)
                ), newX = data.frame(Qn = valid_Q),
                family = stats::gaussian()
              ))
              grn1 <- fm1$pred
              fm2 <- do.call(SL_gr, args = list(
                Y = as.numeric(Aeqa &
                  trainDeltaA == 1 & trainDeltaY == 1),
                X = data.frame(Qn = train_Q), obsWeights = rep(
                  1,
                  length(trainA)
                ), newX = data.frame(Qn = valid_Q),
                family = stats::binomial()
              ))
              grn2 <- fm2$pred
              grn2[grn2 < tolg] <- tolg
            } else if (reduction == "bivariate") {
              fm2 <- do.call(SL_gr, args = list(
                Y = as.numeric(Aeqa &
                  trainDeltaA == 1 & trainDeltaY == 1),
                X = data.frame(Qn = train_Q, gn = train_g),
                obsWeights = rep(1, length(trainA)),
                newX = data.frame(Qn = valid_Q, gn = valid_g),
                family = stats::binomial()
              ))
              grn2 <- fm2$pred
              grn2[grn2 < tolg] <- tolg
              fm1 <- NULL
              grn1 <- rep(NA, length(validY))
            }
          }
        }
        out <- list(grn1 = grn1, grn2 = grn2, fm1 = NULL, fm2 = NULL)
        if (returnModels) {
          out$fm1 <- fm1
          out$fm2 <- fm2
        }
        return(out)
      }
    )
  }

  # GLM
  if (!is.null(glm_gr)) {
    grn <- mapply(
      a = a_0, train_Q = train_Qn, train_g = train_gn,
      valid_Q = valid_Qn, valid_g = valid_gn, SIMPLIFY = FALSE,
      FUN = function(a, train_Q, train_g, valid_Q, valid_g) {
        Aeqa <- trainA == a
        Aeqa[is.na(Aeqa)] <- FALSE

        if (length(unique(train_Q)) == 1) {
          glm_gr <- "1"
        }
        if (reduction == "univariate") {
          fm1 <- stats::glm(
            stats::as.formula(paste0("grn1~", glm_gr)),
            family = "gaussian", data = data.frame(grn1 = (as.numeric(
              Aeqa & trainDeltaA == 1 & trainDeltaY == 1
            ) - train_g) /
              train_g, Qn = train_Q)
          )
          grn1 <- stats::predict(fm1, newdata = data.frame(grn1 = rep(
            0,
            length(validA)
          ), Qn = valid_Q), type = "response")
          fm2 <- stats::glm(
            stats::as.formula(paste0("A~", glm_gr)),
            family = "binomial", data = data.frame(A = as.numeric(Aeqa &
              trainDeltaY == 1 & trainDeltaA == 1), Qn = train_Q)
          )
          grn2 <- stats::predict(
            fm2, type = "response", newdata =
              data.frame(A = rep(0, length(validA)), Qn = valid_Q)
          )
        } else if (reduction == "bivariate") {
          fm1 <- NULL
          grn1 <- rep(NA, length(validY))
          fm2 <- stats::glm(
            stats::as.formula(paste0("A~", glm_gr)),
            family = "binomial", data = data.frame(
              A = as.numeric(Aeqa &
                trainDeltaY == 1 & trainDeltaA == 1), Qn = train_Q,
              gn = train_g
            )
          )
          grn2 <- stats::predict(
            fm2, type = "response",
            newdata = data.frame(
              A = rep(0, length(validA)),
              Qn = valid_Q, gn = valid_g
            )
          )
        }
        grn2[grn2 < tolg] <- tolg
        out <- list(grn1 = grn1, grn2 = grn2, fm1 = NULL, fm2 = NULL)
        if (returnModels) {
          out$fm1 <- fm1
          out$fm2 <- fm2
        }
        return(out)
      }
    )
  }
  tmp1 <- lapply(grn, function(x) {
    data.frame(grn1 = x$grn1, grn2 = x$grn2)
  })
  tmp2 <- lapply(grn, function(x) {
    list(fm1 = x$fm1, fm2 = x$fm2)
  })
  return(list(est = tmp1, fm = tmp2))
}
