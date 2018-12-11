##' The function to initiate the H-family matrix for initial values of
##' x. It includes the x values, and h(x) and h'(x). 
##' @title initialize_hfamily()
##' @param M: Number of observations sampled 
##' @param lb: Lower bound of the domain 
##' @param ub: Upper bound of the domain
##' @param h: Function to compute log-density function
##' @param hprime: Function to compute derivative of log-density function 
##' @param x_initial: Initial values for the x vector 
##' @return A matrix with H-family for initial two values of the x-vector 
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
initialize_hfamily <- function(M, lb, ub, h, hprime, x_initial){
  
  ## Initialize  matrix for h family: the columns of matrix are
  ## x, h(x), h'(x)
  ## Filling in the initial values in matrix
    hfamily <- rbind(c(x_initial[1],
                       h(x_initial[1]),
                       hprime(x = x_initial[1],
                              lb = lb,
                              ub = ub)),
                     c(x_initial[2],
                       h(x_initial[2]),
                       hprime(x = x_initial[2],
                              lb = lb,
                              ub = ub)))
    return(hfamily)
}
