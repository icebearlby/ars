##' .. The given function computes the z values using the initial
##' values for x ..
##' @title initialize_z()
##' @param M: Number of samples 
##' @param x1: Initial parameter I
##' @param x2: Initial parameter II 
##' @param hfamily: Matrix with values of initial values of x and
##'                 related h(x) and h'(x)  
##' @return z: The vector with first three values populated on basis
##'            of x1 and x2 
##' @author Sargam Jain
initialize_z <- function(M, x1, x2, hfamily){

    ## Initialize empty vector for z-values
    z <- rep(NA, as.integer((N^(1/3) + 3)))

    ## Initializing a very small integer limiting to 0
    epsilon = 1e-8
    z[1] <- x1
    z[3] <- x2
    if(z[1] - z[2] < epsilon){
        z[2] <- (hfamily[1, 1] + hfamily[2, 1])/2
    } else {
        z[2] <- ((hfamily[2, 2]-hfamily[1, 2]) -
            (hfamily[2, 1]*hfamily[2, 3]-hfamily[1, 1]*hfamily[1, 3]))/
            (hfamily[1, 3]-hfamily[2, 3])
    }
    return(z)
}
 
