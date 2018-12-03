##' ..The following function is the function to initialize the samples
##' for adaptive rejection sampling given the initial values may/may
##' not provided by the user. The function returns the matrix with
##' h(x), h'(x), and z-vector for the initial parameter values..
##' @title initialize_sample()
##' @param M: Number of samples 
##' @param x1: Initial parameter I
##' @param x2: Initial parameter II 
##' @param width: Band widths of the step taken away/towards estimated
##'               modal value  
##' @param mod_val: The data point at which the density is the highest 
##' @param h: The log-density function 
##' @return ??
##' @author Sargam Jain
initialize_sample <- function(M, x1, x2, width, mod_val, h){

    ## Initializing counter index
    count = 0

    ## When both the initial parameter values are initialized by the user
    if(x1 != -Inf && x2 != Inf){
        hfamily <- initialize_hfamily(M, x1, x2, h, hprime)
        zvalues <- initialize_z(M, x11, x2, hfamily)
    }

    ## When only the upper initial value, x2 is defined by the user
    if(x1 == -Inf && x2 != Inf){
        ## If the mode of the density is higher than the upper bound
        if(mod_val > x2){
            mod_val <- x2 - width
        }
        ## Defining the lower initial value
        x11 <- mod_val
        ## Check if h'(x) exists for the newly defined lower bound
        hprime_check <- compute_deriv(x = x11, func = h)
        while(-Inf < hprime_check && hprime_check <= 0 && count <= 100){
            x11 <- x11 - width
            hprime_check <- compute_deriv(x = x11, func = h)
            count = count + 1
        }
        hfamily <- initialize_hfamily(M, x11, x2, h, hprime)
        zvalues <- initialize_z(M, x1, x22, hfamily)
    }

    ## When only the lower initial value, x1 is defined by the user
    if(x1 != -Inf && x2 == Inf){
        if(mod_val < x1){
            mod_val <- x1 + width
        }
        ## Defining the upper initial value
        x22 <- mod_val
        ## Check if h'(x) exists for the newly defined upper bound
        hprime_check <- compute_deriv(x = x22, func = h)
        while(0 <= hprime_check && hprime_check < Inf && count <= 100){
            x22 <- x22 + width
            hprime_check <- compute_deriv(x = x22, func = h)
            count = count + 1
        }
        hfamily <- initialize_hfamily(M, x1, x22, h, hprime)
        zvalues <- initialize_z(M, x1, x2, hfamily)
    }

    ## When both the initial values are not defined by the user
    if(x1 == -Inf && b == Inf){
        x11 <- mod_val - width
        x22 <- mod_val + width
        ## Check if h'(x) exists for the newly defined bounds
        hprime_check1 <- compute_deriv(x = x11, func = h)
        hprime_check2 <- compute_deriv(x = x22, func = h)
        while(-Inf < hprime_check1 && hprime_check1 <= 0 && count <= 100){
            x11 <- x11 - width
            hprime_check1 <- compute_deriv(x = x11, func = h)
            count = count + 1
        }
        while(0 <= hprime_check2 && hprime_check2 < Inf && count <= 100){
            x22 <- x22 + width
            hprime_check2 <- compute_deriv(x = x22, func = h)
            count = count + 1
        }
        hfamily <- initialize_hfamily(M, x11, x22, h, hprime)
        zvalues <- initialize_z(M, x11, x22, hfamily)

        if(count >= 100) {
            stop("mod_val is invalid, initial point cannot be found. Try another mod_val",
                 .call = FALSE)
        }
    }
    return(list(hfamily, zvalues))
}
