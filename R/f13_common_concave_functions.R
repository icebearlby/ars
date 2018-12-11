## The script holds common log concave functions.

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
dnon_diff<-function(x){
    return(abs(x-3))
}

## c. Non-continous
dnon_conti<-function(x){
    return(dbinom(x, 20, 0.4, log = FALSE))
}
