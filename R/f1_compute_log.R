##' .. The function scales the given density function along the
##' logarithmic scale..
##' @title compute_log 
##' @param x: The data for which distribution is evaluated
##' @param func: The density function for with logarithm needs to be computed
##' @return log_val: The logarithmic value of density function
##' @author Sargam Jain
compute_log <- function(x, func){
    log_val <- log(func(x))
    return(log_val)
}
