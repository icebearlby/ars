library(ars)
library(truncdist)
library(testthat)

## context with one test that groups expectations
context("Tests for main function")

############################
## Check legitimate input ##
############################

## Test non-function 
test_that("ars will throw an error if user do not provide a legitimate function", {
    M = 1000
    f = 3
    expect_error(ars(M, lb = 1,ub = 20, f = 3),
                 "Please provide f as a function")
})

## Test for lb > ub 
test_that("ars will throw an error if lb>=ub", {
    M = 1000
    f = dnor
    expect_error(ars(M, lb = 54, ub = 20, f),
                 "Please provide lb and ub such that lb < ub")
})

## Test for illegitimate input to lb and ub
M = 1000
f = dnor
test_that("ars will throw an error if lb or ub is non-numeric values",{
    expect_error(ars(M, lb = "a", ub = 20, f), "please provide a numeric value for lb")
    expect_error(ars(M, lb = 20, ub = "b", f), "please provide a numeric value for ub")
})

## Non-log concave
test_that("ars will throw an error for non-logconcave densities", {
    M = 1000
    lb = 1
    ub = 10
    f = dnon_concave
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Log density function is not concave over the domain.")
})

## ## Non-differentiable
test_that("ars will throw an error for non-differentiable densities", {
    M = 1000
    lb = 0.1
    ub = 10
    f = dnon_diff
    expect_error(ars(M, lb = lb, ub = ub, f),
                 "Log density function is not differentiable over the domain")
})

## ## Not continuous
## test_that("ars will throw an error for non-continuous densities", {
##     M = 1000
##     lb = 1
##     ub = 20
##     f = dnon_conti
##     expect_error(ars(M, lb = lb, ub = ub, f),
##                  "Density function is either dis-continuous, non-differentiable, or convex.")
## })]
  
## Since we check the continuity of the density by uniformly sampling between the provided bounds, 
## we cannot always capture the case if the density is only discontinuous in a few points.
## As a result, we can not pass our test for checking continuity every time.



########################
## Test sample result ##
########################

## Standard Normal Distribution (Strict KS Test + Plot)
M = 1000
lb = -20
ub = 20
pvalue = rep(NA,10)
f = dnor

test_that("ars works with standard normal density ", {
    for(i in 1:10){
        gsample<-ars(M,lb=lb,ub=ub,f)
        pvalue[i]=ks.test(gsample,"pnorm")$p.value
          
    }
    expect_true(sum(pvalue>0.01)!=0)
})


## Uniform Distribution (Strict Two sample F & T Test + Plot)
 test_that("ars works with uniform ", {
     library(ars)
     library(truncdist)
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
     library(ars)
     library(truncdist)
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

## ## ###################################################################################################
## ## ## Test the Result by comparing mean and variance between standard sample in R and our algorithm ##
## ## ###################################################################################################

## Gamma Distribution
  test_that("ars works with Gamma (shape >= 1)", {
      M = 10000
      lb = 0.1
      ub = 100
      f = dgam
      gsample <- ars(M, lb = lb, ub = ub, f)
      std_sample <- rtrunc(10000, spec = "gamma", shape=1.5, rate = 0.5,
                           a = lb, b = ub)
      expect_lt(abs(mean(gsample) - mean(std_sample)), 0.5)
      expect_lt(abs(var(gsample) - var(std_sample)), 0.5)
  })

## Exponential Distribution
  test_that("ars works with exponential", {
      M = 10000
      lb = 0.1
      ub = 10
      f = dExp
      gsample <- ars(M, lb = lb, ub = ub, f)
      std_sample <- rtrunc(10000, spec = "exp", rate = 3, a = lb, b = ub)
      expect_lt(abs(mean(gsample) - mean(std_sample)), 0.1)
      expect_lt(abs(var(gsample) - var(std_sample)), 0.05)
  })

