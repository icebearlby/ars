##' The function to compute derivative of the density function.
##' @title compute_deriv() 
##' @param x: Domain for the density function 
##' @param func: Log-density function 
##' @param lb: Lower bound of the domain 
##' @param ub: Upper bound of the domain
##' @return A vector with the derivative of log-density function for
##'     the x values provided. 
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
compute_deriv <- function(x, func = compute_log, lb, ub){
    epsilon <- 1e-8
    deriv <- rep(NA, length(x))
    if(sum(x == lb) != 0){
        x_lb <- x[x == lb]
        deriv[x == lb] <- (func(x_lb + epsilon) - func(x_lb)) / epsilon 
    }
    if(sum(x == ub) != 0){
        x_ub <- x[x == ub]
        deriv[x == ub] <- (func(x_ub) - func(x_ub - epsilon)) / epsilon 
    }
    if(sum(is.na(deriv)) != 0){
        deriv[is.na(deriv)] <-
            (func(x[is.na(deriv)]+epsilon) - func(x[is.na(deriv)]-epsilon))/(2*epsilon)
    }
    if (sum(is.na(deriv))!=0){
        stop("Please narrow the boundary,or provide a continuous density.", .call = FALSE)
    }
    return(deriv)
}
