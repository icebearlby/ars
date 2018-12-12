##' The function implements the algorithm for adaptive rejection
##' sampling from log-concave density functions. The algorithm
##' implemented is a replication of the methodology proposed in Gilks
##' et al (1992). 
##' 
##' .. content for \details{} ..
##' @title ars() 
##' @param M: Number of observations sampled
##' @param lb: Lower bound for the domain (by default -Inf)
##' @param ub: Upper bound for the domain (by default +Inf)
##' @param f: Function to compute univariate log-concave density function 
##' @param width: Step value around highest probability density point 
##' @param mod_val: Highest probability density point 
##' @return A vector of final M samples generated from adaptive rejection algorithm.
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
ars <- function(M, lb = -Inf, ub = Inf, f, width = 0.5, mod_val = 0){

    ## Check the class of f: Function
    if(class(f) != "function"){
        stop("Please provide f as a function", call. = FALSE)
    }
    #Check if lb and ub are both numeric values
    if(class(lb)!="numeric"){
    stop('please provide a numeric value for lb', call. = FALSE)
    }
    if(class(ub)!="numeric"){
    stop('please provide a numeric value for ub', call. = FALSE)
    }
    ## Check if user provide lb<ub
    if(lb >= ub){
        stop("Please provide lb and ub such that lb < ub", call. = FALSE)
    }
    
    ## Generate samples for sanity check: concavity, differentiability, continuity
    if (lb != -Inf && ub != Inf){
        check_vector <- runif(100, lb, ub)
    }
    if (lb == -Inf && ub != Inf){
        neg_epsilon = -(1e+4)
        check_vector <- runif(100, neg_epsilon, ub)
    }
    if (lb != -Inf && ub == Inf){
        pos_epsilon = 1e+4
        check_vector <- runif(100, lb, pos_epsilon)
    }
    if (lb == -Inf && ub == Inf){
        pos_epsilon = 1e+4
        neg_epsilon = -(1e+4)
        check_vector <- runif(100, neg_epsilon, pos_epsilon)
    }
    check_vector <- sort(check_vector)
    epsilon = 1e-6

    ## Condition for continuity
    if(sum(round(compute_log(check_vector + epsilon), digits = 6) !=
           round(compute_log(check_vector + epsilon), digits = 6)) != 0){
        stop("Density function is either dis-continuous, non-differentiable, or convex.",
    .call = FALSE)
    }
    
    ## Condition for differentiability
    if(sum(round(compute_deriv(x = check_vector,
                               lb = lb,
                               ub = ub),
                 digits=4) != round(compute_deriv(x = (check_vector +
                                                       1e-16),
                                                  lb = lb,
                                                  ub = ub),
                                    digits=4)) != 0){
        stop("Density function is either dis-continuous, non-differentiable, or convex.",
             .call = FALSE)
    }
    
    ## Condition for concavity
    if(sum(round(diff(compute_deriv(x = check_vector, lb = lb, ub = ub)),digits=8) > 1e-4) !=0){
        stop("Density function is either dis-continuous, non-differentiable, or convex.",
             .call = FALSE)
    }
    
    ## Special case for Uniform
    if(sum(compute_deriv(x = check_vector,lb = lb,ub = ub) != 0) == 0){
        gsample = runif(M, min = lb, max = ub)
        hist(gsample)
        return(gsample)
    }

    ## StageI: Initialize sample
    initialized_sample <- initialize_sample(M, lb, ub, width = 0.5, mod_val = 0)
    hfamily <- initialized_sample$hfamily
    zvalues <- initialized_sample$zvalues
    
    sample_count = 0
    sample_bag = 0
    
    while(sample_count <= M){

        ## Stage II: Sampling
        samples <- create_samples(M1 = M^(2/3)+10, hfamily, zvalues)
        ## since there would be around M^(1/3), on average
        ## each time we need the squeezing test to accept M^(2/3) samples

        ## Stage III: Squeezing and Rejection
        sr_test_x <- sr_test(samples, hfamily, zvalues)
        sample_count = sample_count + sr_test_x$count_accept
        sample_bag = c(sample_bag, sr_test_x$x_accept)

        ## Stage IV: Updation
        if (is.na(sr_test_x$x_r) != TRUE){
            hfamily <- update_hfamily(hfamily,x_r = sr_test_x$x_r)
            zvalues <- update_z(hfamily,lb,ub)
        }
        
    }

    ## Final check for concavity
    if(sum(round(diff(hfamily[,3]),digits=8) > epsilon) !=0){
        stop("Density function is either dis-continuous, non-differentiable, or convex.",
             .call = FALSE)
    }  

    ## Return the M sample in the sample bag
    hist(sample_bag[1:M])
    return(sample_bag[1:M])
}
