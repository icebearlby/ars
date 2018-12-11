##' Function to update the H-family matrix with new samples i.e. the
##' first sample that did not pass the squeezing test.
##' @title update_hfamily()
##' @param hfamily: H-family matrix 
##' @param x_r: First sample that did not pass the squeezing test 
##' @param h: Function to compute log-density function 
##' @param hprime: Function to compute derivative of log-density function 
##' @return A matrix with the updated H-Family
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
update_hfamily <- function(hfamily, x_r, h = compute_log, hprime = compute_deriv){
  
    newline = cbind(x_r, h(x_r), hprime(x_r, lb = lb, ub = ub))
    hfamily = rbind(hfamily, newline)
    hfamily = hfamily[order(hfamily[ , 1]), ]
    return(hfamily)
}
