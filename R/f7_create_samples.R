##' The function to generate samples by converting uniform samples to
##' X-vector according to sampled CDF.
##' @title create_samples()
##' @param M1: Number of uniform samples generated in each iteration 
##' @param hfamily: H-family matrix 
##' @param zvalues: Z-vector 
##' @param scdf: Function to compute sampled CDF  
##' @return A vector of samples generated
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
create_samples <- function(M1, hfamily, zvalues, scdf = calculate_scdf){
    zcdf <- calculate_scdf(vals = zvalues, hfamily, z = zvalues)
    zq = zcdf$scdf
    c = zcdf$c
    unif_sample = runif(M1)
    uidx = findInterval(unif_sample, zq)
    intervals_count = length(zq) - 1
    zlow = zvalues[-length(zvalues)]
    res = rep(NA, length(unif_sample))
    
    for(ii in 1:intervals_count){
        ui = unif_sample[uidx == ii]
        if(length(ui) == 0){next}
        
        xp = hfamily[(ii),1]
        hp = hfamily[(ii),2]
        hprimep = hfamily[(ii),3]
        zl = zlow[(ii)]
        ## invert the CDF
        tmp = log((ui-zq[ii]) * hprimep * c / exp(hp) + exp((zl - xp)*hprimep))/hprimep + xp
        ## convert uniform sample to x according to sampled cdf
        res[uidx == ii] = tmp
    }  
    return(res)
}
