##' .. The following function initializes the values of h(x) and h'(x)
##' given the initial two parameter values for x ..
##' @title initialize_hfamily()
##' @param M: Number of samples
##' @param x1: Initial parameter I 
##' @param x2: Initial parameter II
##' @param h: Log density function 
##' @param hprime: Derivative of log differenced function
##' @return hfamily: Matrix with x-values, h(x), and h'(x) defined for
##'                  initial parameter values 
##' @author Sargam Jain
initialize_hfamily <- function(M, x1, x2, h, hprime){

    ## Initialize empty matrix for h family: the columns of matrix are
    ## x, h(x), h'(x)
    hfamily <- matrix(NA,
                      nrow = as.integer((M^(1/3) + 2)),
                      ncol = 3)
    
    ## Filling in the initial values in matrix
    hfamily[1, 1] <- x1
    hfamily[1, 2] <- h(x1)
    hfamily[1, 3] <- hprime(x1)
    hfamily[2, 1] <- x2
    hfamily[2, 2] <- h(x2)
    hfamily[2, 3] <- hprime(x2)

    ## Tests for finite values for log density function
    assert_that(hfamily[1, 2])
    assert_that(hfamily[2, 2])

    ## Test for finite values for derivative of log density function
    assert_that(hfamily[1, 3])
    assert_that(hfamily[2, 3])

    return(hfamily)
}

