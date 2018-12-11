##' The function to perform the squeezing test and rejection test. 
##' @title sr_test() 
##' @param x_sample: Samples generated in the sampling stage  
##' @param hfamily: H-Family matrix 
##' @param zvalues: Z-vector 
##' @param h: Function to compute log-density function 
##' @return A list with samples accpeted, first sample that did not
##'     pass the squeezing test, and counter for number of samples accepted.
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
sr_test <- function(x_sample, hfamily, zvalues, h = compute_log){
  
    w = runif(length(x_sample))
    count_accept = 0
    x_r = NA
  
    ## perform squeezing test and select the samples that passed the squeezing test
    x_accept = x_sample[exp(l(x_sample,hfamily) - u(x_sample,hfamily,zvalues)) >= w]
    x_reject_s = x_sample[exp(l(x_sample,hfamily) - u(x_sample,hfamily,zvalues)) < w]
    count_accept = count_accept + length(x_accept)
    
    if(length(x_reject_s) != 0){
        ## perform rejection test and select the first sample that did not passed the squeezing test
        x_r = x_reject_s[1]
        
        if(exp(h(x_reject_s[1]) - u(x_reject_s[1],hfamily,zvalues)) >= w[x_sample==x_r]){
            count_accept = count_accept + 1
            x_accept = c(x_accept,x_r)
        }
    }
    list = list(x_accept = x_accept, x_r = x_r, count_accept = count_accept)
    return(list)
}
