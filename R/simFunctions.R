#' makeDataBiv
#' 
#' Make simulation data set of size n
#' @param n Sample size
#' @param beta Strength of interaction
#' 
#' @importFrom stats plogis rbinom runif

makeDataBiv <- function(n, beta = 1){
    W <- data.frame(W1=stats::runif(n,-2,2), W2=stats::rbinom(n, 1, 0.5))
    A <- stats::rbinom(n, 1, stats::plogis(-beta*W$W1 + 2*beta*W$W1*W$W2))
    Y <- stats::rbinom(n, 1, stats::plogis(0.2*A - beta*W$W1 + 2*beta*W$W1*W$W2))

    return(list(W=W, A=A, Y=Y))
}

#' makeDataSix
#' Make simulation data set of size n
#' @param n Sample size
#' @param beta Strength of interaction
#' @importFrom stats plogis rbinom runif
#' 

makeDataSix <- function(n, beta = 1){
    W <- data.frame(W1=stats::runif(n,-2,2), W2=stats::rbinom(n, 1, 0.5), W3 = stats::runif(n,-2,2), W4 = stats::rbinom(n,1,0.5), 
                    W5=stats::runif(n,-2,2), W6 = stats::rbinom(n,1,0.5))
    A <- stats::rbinom(n, 1, stats::plogis(-beta*W$W1 + 2*beta*W$W1*W$W2))
    Y <- stats::rbinom(n, 1, stats::plogis(0.2*A - beta*W$W1 + 2*beta*W$W1*W$W2))

    return(list(W=W, A=A, Y=Y))
}


#' getTruth
#' get true parameter value for simulation
#' @param beta Strength of interaction 
#' @param a0 Set level of treatment
#' 
getTruth <- function(beta=1,a0=1){
  # here the -2/2 comes from the limits of the uniform distribution
  # this is from wolfram, has been checked numerically for accuracy
  1/2 + 1/(8*beta) * log((exp(0.2*a0) + exp(-2*beta))/(1+exp(0.2*a0 - 2*beta)) /((exp(0.2*a0) + exp(2*beta))/(1+exp(0.2*a0 + 2*beta))))
}

