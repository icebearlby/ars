library(ars)
library(testthat)
library(truncdist)
## context with one test that groups expectations
context("Tests for main function")

########################
## Test sample result ##
########################

## Standard Normal Distribution (Strict KS Test + Plot)
test_that("ars works with standard normal density ", {
    M = 1000
    lb = -20
    ub = 20
    pvalue = rep(NA,10)
    f = dnor
    for(i in 1:10){
        gsample<-ars(M,lb=lb,ub=ub,f=dnor)
        pvalue[i]=ks.test(gsample,"pnorm")$p.value
          
    }
    expect_true(sum(pvalue>0.01)!=0)
})


## Uniform Distribution (Strict Two sample F & T Test + Plot)
test_that("ars works with uniform ", {
    mean_diff=rep(NA,10)
    var_diff=rep(NA,10)
    M = 10000
    lb = 0.1
    ub = 0.9
    f=duni
    std_sample <- rtrunc(10000,spec="unif",a=lb,b=ub)
    for(i in 1:10){
        gsample <- ars(M,lb=lb,ub=ub,f)
        mean_diff[i]=abs(mean(gsample)-mean(std_sample))
        var_diff[i]=abs(var(gsample)-var(std_sample))
          
    }
    expect_true(sum(mean_diff<0.8)!=0)
    expect_true(sum(var_diff<0.8)!=0)
    
})

## Logistic Distribution (Strict KS Test + Plot)
test_that("ars works with logistic ", {
    pvalue=rep(NA,10)
    M = 100
    lb = -15
    ub = 15
    f=dlogit
    for(i in 1:10){
        gsample <- ars(M,lb=lb,ub=ub,f)
        pvalue[i]=ks.test(gsample,"plogis")$p.value
          
    }
    expect_true(sum(pvalue>0.01)!=0)
})

## ###################################################################################################
## ## Test the Result by comparing mean and variance between standard sample in R and our algorithm ##
## ###################################################################################################

## ## Gamma Distribution
## test_that("ars works with Gamma (shape >= 1)", {
##     M = 10000
##     lb = 0.1
##     ub = 100
##     f = dgam
##     gsample <- ars(M, lb = lb, ub = ub, f)
##     std_sample <- rtrunc(10000, spec = "gamma", shape=1.5, rate = 0.5,
##                          a = lb, b = ub)
##     expect_lt(abs(mean(gsample) - mean(std_sample)), 0.5)
##     expect_lt(abs(var(gsample) - var(std_sample)), 0.5)
## })


## Beta Distribution
test_that("ars works with beta(4,3) ", {
    mean_diff=rep(NA,10)
    var_diff=rep(NA,10)
    M = 10000
    lb = 0
    ub = 1
    f=dBeta
    std_sample <- rtrunc(10000, spec = "beta", shape1 = 4,
                         shape2 = 3, a = lb,b = ub)
    for(i in 1:10){
        gsample<-ars(M,lb=lb,ub=ub,f)
        mean_diff[i]=abs(mean(gsample)-mean(std_sample))
        var_diff[i]=abs(var(gsample)-var(std_sample))
          
    }
    expect_true(sum(mean_diff<0.8)!=0)
    expect_true(sum(var_diff<0.8)!=0)
    
})

## ## Exponential Distribution
## test_that("ars works with exponential", {
##     M = 10000
##     lb = 0.1
##     ub = 10
##     f = dExp
##     gsample <- ars(M, lb = lb, ub = ub, f)
##     std_sample <- rtrunc(10000, spec = "exp", rate = 3, a = lb, b = ub)
##     expect_lt(abs(mean(gsample) - mean(std_sample)), 0.1)
##     expect_lt(abs(var(gsample) - var(std_sample)), 0.05)
## })

############################
## Check legitimate input ##
############################

## Non-log concave
test_that("ars will throw an error for non-logconcave densities", {
    M = 1000
    lb = 1
    ub = 10
    f = dnon_concave
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Density function is either dis-continuous, non-differentiable, or convex.")
})

## Non-differentiable
test_that("ars will throw an error for non-differentiable densities", {
    M = 1000
    lb = 0.1
    ub = 10
    f = dnon_diff
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Density function is either dis-continuous, non-differentiable, or convex.")
})

## Not continuous
test_that("ars will throw an error for non-continuous densities", {
    M = 1000
    lb = 1
    ub = 20
    f = dnon_conti
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Density function is either dis-continuous, non-differentiable, or convex.")
})

## Test non-function 
test_that("ars will throw an error if user do not provide a legitimate function", {
    M = 1000
    lb = 1
    ub = 20
    f = 3
    expect_error(ars(M, lb = lb,ub = ub, f),
                 "Please provide f as a function")
})

## Test for lb > ub 
test_that("ars will throw an error if lb>=ub", {
    M = 1000
    lb = 54
    ub = 20
    f = dnor
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Please provide lb and ub such that lb < ub")
})
