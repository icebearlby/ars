##' The function to initialize the H-family and Z-vector under four
##' scenarios, designed on basis of whether the bounds have been
##' provided by the user.
##' @title initialize_sample()
##' @param M: Number of observations sampled 
##' @param lb: Lower bound of the domain
##' @param ub: Upper bound of the domain 
##' @param width: Step value around highest probability density point 
##' @param mod_val: Highest probability density point  
##' @param h: Function to compute log-density function 
##' @param hprime: Function to compute derivative of log-density function  
##' @return A list of H-family matrix and Z-vector
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
initialize_sample <- function(M, lb, ub, width = 0.5, mod_val = 0,
                              h = compute_log,
                              hprime = compute_deriv){

    ## Initializing counter index
    count = 0

    ## When both the initial parameter values are initialized by the user
    if(lb != -Inf && ub != Inf){
        x_initial <- sort(runif(2, lb, ub))
        hfamily <- initialize_hfamily(M, lb, ub, h, hprime, x_initial)
        zvalues <- initialize_z(M, lb, ub, hfamily)
    }

    ## When only the upper initial value, x2 is defined by the user
    if(lb == -Inf && ub != Inf){
        ## If the mode of the density is higher than the upper bound
        if(mod_val > ub){
            mod_val <- ub - width
        }
        ## Defining the lower initial value
        llb <- mod_val
        ## Check if h'(x) exists for the newly defined lower bound
        hprime_check <- hprime(x = llb, func = h, lb = lb, ub = ub)
        while(-Inf < hprime_check && hprime_check <= 0 && count <= 100){
            llb <- llb - width
            hprime_check <- compute_deriv(x = llb, func = h, lb = lb, ub = ub)
            count = count + 1
        }
        x_initial <- sort(c(llb, runif(1, llb, ub)))
        hfamily <- initialize_hfamily(M, llb, ub, h, hprime,x_initial)
        zvalues <- initialize_z(M, lb, ub, hfamily)
    }

    ## When only the lower initial value, x1 is defined by the user
    if(lb != -Inf && ub == Inf){
        if(mod_val < lb){
            mod_val <- lb + width
        }
        ## Defining the upper initial value
        uub <- mod_val
        ## Check if h'(x) exists for the newly defined upper bound
        hprime_check <- compute_deriv(x = uub, func = h, lb = lb, ub = ub)
        while(0 <= hprime_check && hprime_check < Inf && count <= 100){
            uub <- uub + width
            hprime_check <- compute_deriv(x = uub, func = h, lb = lb, ub = ub)
            count = count + 1
        }
        x_initial <- sort(c(runif(1, lb, uub),uub))
        hfamily <- initialize_hfamily(M, lb, uub, h, hprime,x_initial)
        zvalues <- initialize_z(M, lb, ub, hfamily)
    }

    ## When both the initial values are not defined by the user
    if(lb == -Inf && ub == Inf){
        llb <- mod_val - width
        uub <- mod_val + width
        ## Check if h'(x) exists for the newly defined bounds
        hprime_check1 <- compute_deriv(x = llb, func = h, lb = lb, ub = ub)
        hprime_check2 <- compute_deriv(x = uub, func = h, lb = lb, ub = ub)
        while(-Inf < hprime_check1 && hprime_check1 <= 0 && count <= 100){
            llb <- llb - width
            hprime_check1 <- compute_deriv(x = llb, func = h, lb = lb, ub = ub)
            count = count + 1
        }
        while(0 <= hprime_check2 && hprime_check2 < Inf && count <= 100){
            uub <- uub + width
            hprime_check2 <- compute_deriv(x = uub, func = h, lb = lb, ub = ub)
            count = count + 1
        }
        x_initial <- c(llb,uub)
        hfamily <- initialize_hfamily(M, llb, uub, h, hprime,x_initial)
        zvalues <- initialize_z(M, lb, ub, hfamily)

        if(count >= 100) {
            stop("mod_val is invalid, initial point cannot be found. Try another mod_val",
                 .call = FALSE)
        }
    }
    return(list("hfamily" = hfamily, "zvalues" = zvalues))
}
