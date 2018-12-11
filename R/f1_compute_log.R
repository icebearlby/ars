##' Function to compute log-density function from the density function.
##' @title compute_log()
##' @param x: Domain for density function 
##' @param func: Density function 
##' @return A vector with values of log density function
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
compute_log <- function(x, func = f){
    log_val <- log(func(x))
    return(log_val)
}
