##' Function to compute log-density function from the density function.
##' @title compute_log()
##' @param x: Domain for density function 
##' @param func: Density function 
##' @return A vector with values of log density function
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
compute_log <- function(x, func = f){
    log_val <- log(func(x))
    return(log_val)
}


##' The function to compute derivative of the density function.
##' @title compute_deriv() 
##' @param x: Domain for the density function 
##' @param func: Log-density function 
##' @param lb: Lower bound of the domain 
##' @param ub: Upper bound of the domain
##' @return A vector with the derivative of log-density function for
##'     the x values provided. 
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
compute_deriv <- function(x, func = compute_log, lb, ub, f){
    epsilon <- 1e-8
    deriv <- rep(NA, length(x))
    if(sum(x == lb) != 0){
        x_lb <- x[x == lb]
        deriv[x == lb] <- (func(x_lb + epsilon) - func(x_lb)) / epsilon 
    }
    if(sum(x == ub) != 0){
        x_ub <- x[x == ub]
        deriv[x == ub] <- (func(x_ub) - func(x_ub - epsilon)) / epsilon 
    }
    if(sum(is.na(deriv)) != 0){
        deriv[is.na(deriv)] <- (func(x[is.na(deriv)]+epsilon)
            - func(x[is.na(deriv)]))/ epsilon
    }
    if (sum(is.na(deriv))!=0){
        stop("Please narrow the boundary, or provide a continous density.",
             .call = FALSE)
    }
    return(deriv)
}


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


##' The function to compute the total CDF and CDF at each interval.
##' @title calculate_scdf
##' @param val: Z-vector
##' @param hfamily: H-family matrix 
##' @param z: Z-vector 
##' @return A list with CDF at each interval and total CDF
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
calculate_scdf <- function(vals, hfamily, z) {
  
    zlen = length(z)
    cdf = numeric(length(vals))
    c = 0
    
    for(i in 1:(zlen-1)){
        zl = z[i]
        zu = z[i+1]
        xp = hfamily[i,1]
        hp = hfamily[i,2]
        hprimep = hfamily[i,3]
        ## calculate the cumulated density in each interval, parts of the demonimator of s
        ds = exp(hp)/hprimep * ( exp((zu - xp)*hprimep) - exp((zl - xp)*hprimep) )
        ## get the logic vector, true if val belongs to the
        ## interval
        inside_idx = (zl < vals & vals <= zu)
        ## get the logic vector, true if val is larger than z[i+1]
        greater_idx = vals > zu
        ## if inside interval, only cumulate till the example
        ## point
        cdf[inside_idx] = cdf[inside_idx] + exp(hp)/hprimep * (exp((vals[inside_idx] - xp)*hprimep) - exp((zl - xp)*hprimep))
        ## if larger than the upper interval, cumulate the whole
        ## interval
        cdf[greater_idx] = cdf[greater_idx] + ds
        ## total CDF
        c = c + ds
    }
    l = list(scdf = cdf/c, c = c )
    return(l)
}


##' The function to generate samples by converting uniform samples to
##' X-vector according to sampled CDF.
##' @title create_samples()
##' @param M1: Number of uniform samples generated in each iteration 
##' @param hfamily: H-family matrix 
##' @param zvalues: Z-vector 
##' @param scdf: Function to compute sampled CDF  
##' @return A vector of samples generated
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
create_samples <- function(M1, hfamily, zvalues, scdf = calculate_scdf){
    zcdf <- calculate_scdf(vals = zvalues, hfamily, z = zvalues)
    zq = zcdf$scdf
    c = zcdf$c
    unif_sample = runif(M1)
    uidx = findInterval(unif_sample, zq)
    intervals_count = length(zq) - 1
    zlow = zvalues[-length(zvalues)]
    res = rep(NA, length(unif_sample))
    
    for(ii in 1:intervals_count){
        ui = unif_sample[uidx == ii]
        if(length(ui) == 0){next}
        
        xp = hfamily[(ii),1]
        hp = hfamily[(ii),2]
        hprimep = hfamily[(ii),3]
        zl = zlow[(ii)]
        ## invert the CDF
        tmp = log((ui-zq[ii]) * hprimep * c / exp(hp) + exp((zl - xp)*hprimep))/hprimep + xp
        ## convert uniform sample to x according to sampled cdf
        res[uidx == ii] = tmp
    }  
    return(res)
}


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


##' The function to compute the lower hull used to generate squeezing function.
##' @title l() 
##' @param x_sample: Samples generated in sampling stage 
##' @param hfamily: H-family matrix
##' @return A vector of values for l(x)
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
l <- function(x_sample, hfamily){
    res = rep(0, length(x_sample))
    interval_idx = findInterval(x_sample, hfamily[,1])
    
    nintervals = length(hfamily[,1]) - 1
    for(idx in 1:nintervals){
        xi = x_sample[interval_idx == idx]
        xx  = ( (hfamily[idx+1,1] - xi)*hfamily[idx,2] + (xi - hfamily[idx,1])*hfamily[idx+1,2] ) / ( hfamily[idx+1,1] - hfamily[idx,1])
        res[interval_idx == idx] = xx
    }
    
    res[interval_idx == 0] = -Inf
    res[interval_idx == length(hfamily[,1])]=-Inf
    return(res)
}


##' The function to perform the squeezing test and rejection test. 
##' @title sr_test() 
##' @param x_sample: Samples generated in the sampling stage  
##' @param hfamily: H-Family matrix 
##' @param zvalues: Z-vector 
##' @param h: Function to compute log-density function 
##' @return A list with samples accpeted, first sample that did not
##'     pass the squeezing test, and counter for number of samples accepted.
##' @author Baoyue Liang, Sargam Jain, Dandan Ru
sr_test <- function(x_sample, hfamily, zvalues, h = compute_log){
  
    w = runif(length(x_sample))
    count_accept = 0
    x_r = NA
  
    ## perform squeezing test and select the samples that passed the squeezing test
    x_accept = x_sample[exp(l(x_sample,hfamily) - u(x_sample,hfamily,zvalues)) >= w]
    x_reject_s = x_sample[exp(l(x_sample,hfamily) - u(x_sample,hfamily,zvalues)) < w]
    count_accept = count_accept + length(x_accept)
    
    if(length(x_reject_s) != 0){
        ## perform rejection test and select the first sample that did not passed the squeezing test
        x_r = x_reject_s[1]
        
        if(exp(h(x_reject_s[1]) - u(x_reject_s[1],hfamily,zvalues)) >= w[x_sample==x_r]){
            count_accept = count_accept + 1
            x_accept = c(x_accept,x_r)
        }
    }
    list = list(x_accept = x_accept, x_r = x_r, count_accept = count_accept)
    return(list)
}


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

############################### Non-legit functions #########################################

########################
## 1. Standard normal ##
########################
dnor <- function(x){ 
    return((1/(sqrt(2*pi)))*exp(-(x^2)/2))
}


############################
## 2. Gamma (chi-squared) ##
############################
dgam <- function(x){
    return((1/(gamma(3/2)*(1/2)^(3/2))*(x^(1/2))*exp(-x*1/2)))
}


#####################
## 3. Uniform(0,1) ##
#####################
duni <- function(x){
    return(dunif(x, min=0, max=1))
}


##################
## 4. Beta(4,3) ##
##################
dBeta <- function(x){
    return(dbeta(x, shape1 = 4, shape2 = 3))
}


#################
## 5. Logistic ##
#################
dlogit <- function(x){
    return(dlogis(x))
}

################
## 6. T-Dist  ##
################
d_t <- function(x){
    return(dt(x, df = 1, ncp = 0, log = FALSE))
}


################
## 7. F-Dist  ##
################
dF <- function(x){
    return(df(x, 10, 1, ncp = 0, log = FALSE))
}


#####################
## 8. Exponential  ##
#####################
dExp <- function(x){
    return(dexp(x, rate = 3, log = FALSE))
}


###################
## 9. Lognormal  ##
###################
dlnor <- function(x){
    return(dlnorm(x))
}

#########################
## Non-legit Functions ##
#########################

## a. Non-concave 
dnon_concave <- function(x){
    return(exp(x^2))
}

## b. Non-differentiable
dnon_diff <- function(x){
    return(exp(abs(x-3)))
}

## c. Non-continous
dnon_conti<-function(x){
    return(dbinom(x, 20, 0.4, log = FALSE))
}





