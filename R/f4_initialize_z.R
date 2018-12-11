##' The function to initialize the z-vector i.e the points at which
##' tangents to x points intersect
##' @title initialize_z()
##' @param M: Number of observations sampled
##' @param lb: Lower bound of the domain 
##' @param ub: Upper bound of the domain 
##' @param hfamily: H-family matrix 
##' @return A vector with the initial three values of the z-vector
##' @author Sargam Jain
initialize_z <- function(M, lb, ub, hfamily){
  z<-rep(NA,3)
  ## Initializing a very small integer limiting to 0
  epsilon = 1e-8
  z[1] <- lb
  z[3] <- ub
  if(abs(z[1] - z[3]) < epsilon){
    z[2] <- (hfamily[1, 1] + hfamily[2, 1])/2
  } else {
    z[2] <- ((hfamily[2, 2]-hfamily[1, 2]) -
               (hfamily[2, 1]*hfamily[2, 3]-hfamily[1, 1]*hfamily[1, 3]))/
      (hfamily[1, 3]-hfamily[2, 3])
  }
  if(is.infinite(z[2])){
  z[2] = ( hfamily[1,1] + hfamily[2,1] )/2
  }
  return(z)
}
