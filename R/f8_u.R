##' The function to compute the upper hull used to generate envelope function.
##' @title u()
##' @param x_sample: Samples generated in sampling stage 
##' @param hfamily: H-family matrix 
##' @param z: Z-vector 
##' @return A vector of values for u(x)
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
u <- function(x_sample, hfamily,z) {
  
  res = rep(0, length(x_sample))
  interval_idx = findInterval(x_sample, z)
  
  nintervals = length(z) -1
  for(idx in 1:nintervals){
    xi = x_sample[interval_idx == idx]
    ux = hfamily[idx,2] + (xi - hfamily[idx,1]) * hfamily[idx,3]
    res[interval_idx == idx] = ux
  }
  
  return(res)
}
