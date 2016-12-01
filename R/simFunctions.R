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