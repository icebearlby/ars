##' The function to compute the lower hull used to generate squeezing function.
##' @title l() 
##' @param x_sample: Samples generated in sampling stage 
##' @param hfamily: H-family matrix
##' @return A vector of values for l(x)
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
l <- function(x_sample, hfamily){
    res = rep(0, length(x_sample))
    interval_idx = findInterval(x_sample, hfamily[,1])
    
    nintervals = length(hfamily[,1]) - 1
    for(idx in 1:nintervals){
        xi = x_sample[interval_idx == idx]
        xx  = ( (hfamily[idx+1,1] - xi)*hfamily[idx,2] + (xi - hfamily[idx,1])*hfamily[idx+1,2] ) / ( hfamily[idx+1,1] - hfamily[idx,1])
        res[interval_idx == idx] = xx
    }
    
    res[interval_idx == 0] = -Inf
    res[interval_idx == length(hfamily[,1])]=-Inf
    return(res)
}
