##' .. The following function computes the derivative of the function
##' provided as argument to the function..
##'
##' @title compute_deriv 
##' @param x: data provided for the analysis 
##' @param func: the log-density for which derivative needs to be
##'     computed 
##' @param x1: Initial values provided by user for adaptive rejection sampling
##' @param x2: Initial values provided by user for adaptive rejection sampling
##' @return deriv: the derivative of the log density function 
##' @author Sargam Jain
compute_deriv <- function(x, func, x1, x2){
    epsilon <- 1e-8
    ## Given that x1 and x2 are bounds of the function:
    if(x == x1){
        deriv <- (func(x + epsilon) - func(x)) / epsilon 
    }
    if (x == x2){
        deriv <- (func(x) - func(x - epsilon)) / epsilon
    }
    if(x >= x1 && x <= x2){
        deriv <- (func(x + epsilon) - func(x - epsilon))/(2*epsilon)
    }
    return(deriv)
}
 
