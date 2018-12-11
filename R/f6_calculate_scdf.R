##' The function to compute the total CDF and CDF at each interval.
##' @title calculate_scdf
##' @param val: Z-vector
##' @param hfamily: H-family matrix 
##' @param z: Z-vector 
##' @return A list with CDF at each interval and total CDF
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
calculate_scdf <- function(vals, hfamily, z) {
  
    zlen = length(z)
    cdf = numeric(length(vals))
    c = 0
    
    for(i in 1:(zlen-1)){
        zl = z[i]
        zu = z[i+1]
        xp = hfamily[i,1]
        hp = hfamily[i,2]
        hprimep = hfamily[i,3]
        ## calculate the cumulated density in each interval, parts of the demonimator of s
        ds = exp(hp)/hprimep * ( exp((zu - xp)*hprimep) - exp((zl - xp)*hprimep) )
        ## get the logic vector, true if val belongs to the
        ## interval
        inside_idx = (zl < vals & vals <= zu)
        ## get the logic vector, true if val is larger than z[i+1]
        greater_idx = vals > zu
        ## if inside interval, only cumulate till the example
        ## point
        cdf[inside_idx] = cdf[inside_idx] + exp(hp)/hprimep * (exp((vals[inside_idx] - xp)*hprimep) - exp((zl - xp)*hprimep))
        ## if larger than the upper interval, cumulate the whole
        ## interval
        cdf[greater_idx] = cdf[greater_idx] + ds
        ## total CDF
        c = c + ds
    }
    l = list(scdf = cdf/c, c = c )
    return(l)
}
