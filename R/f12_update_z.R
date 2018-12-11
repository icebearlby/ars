##' Function to update the Z-vector using the updated H-family.
##' @title update_z() 
##' @param hfamily: Updated H-family matrix 
##' @param lb: Lower bound of the domain 
##' @param ub: Upper bound of the domain 
##' @return A vector with updated Z-values
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
update_z <- function(hfamily, lb, ub){
    nrow = nrow(hfamily)
    xf0 = hfamily[-nrow, 1]
    xf1 = hfamily[-1, 1]
    hf0 = hfamily[-nrow, 2]
    hf1 = hfamily[-1, 2]
    dhf0 = hfamily[-nrow, 3]
    dhf1 = hfamily[-1, 3]
    
    z = xf0 + (hf0 - hf1 + (xf1 - xf0)*dhf1) / (dhf1 - dhf0)
    inf_idx = is.infinite(z) == TRUE
    z[inf_idx] = ( xf1[inf_idx] + xf0[inf_idx] )/2
    
    z = c(lb, z, ub)
    z = sort(z, decreasing = FALSE)
    return(z)	
}
