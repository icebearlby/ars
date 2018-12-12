library(ars)
library(testthat)
context("Testing the Utility Functions")

###########################
## Test: compute_deriv() ##
###########################
## lb = -100
## ub = 100
## f = dnor
## test_that("compute_deriv will throw an error if the boundary is too wide for certain density", {
##     library(ars)
##     expect_error(compute_deriv(99, lb = lb,ub = ub),
##                  "Please narrow the boundary, or provide a continous density.")
##     expect_error(compute_deriv(-99,lb = lb, ub = ub),
##                  "Please narrow the boundary, or provide a continous density.")
## })


## ###############################
## ## Test: initialize_sample() ##
## ###############################
## M = 1000
## f = dnor
## test_that("Initialize step choose proper starting points", {
##     library(ars)
##     expect_gt(initialize_sample(M, -Inf, 8)$hfamily[1, 3], 0)
##     expect_lt(initialize_sample(M, -8, Inf)$hfamily[2, 3], 0)
##     expect_gt(initialize_sample(M, -Inf, Inf)$hfamily[1, 3], 0)
##     expect_lt(initialize_sample(M, -Inf, Inf)$hfamily[2, 3], 0)
## })

## ## Test x in the boundary


## ##########################
## ## Test: initialize_z() ##
## ##########################
## M = 1000
## lb = 0.1
## ub = 30
## f = dexp
## test_that("initialize_z do not return Inf values", {
##     library(ars)
##     expect_equal(sum(is.infinite(initialize_sample(M, lb, ub)$zvalues)),
##                  0)
## })

## #################################################
## ## Test: update_z() (do not return Inf values) ##
## #################################################
## M = 1000
## lb = -10
## ub = 10
## f = dnor
## test_that("update_z do not return Inf values", {
##     library(ars)
##     initialized_sample <- initialize_sample(M, lb, ub, width = 0.5, mod_val = 0)
##     hfamily <- initialized_sample$hfamily
##     zvalues <- initialized_sample$zvalues
##     while(TRUE){
##         samples <- create_samples(M1 = M^(2/3)*2, hfamily, zvalues)
##         sr_test_x <- sr_test(samples,hfamily,zvalues)
##         if (is.na(sr_test_x$x_r) != TRUE){
##             hfamily <- update_hfamily(hfamily,x_r = sr_test_x$x_r)
##             zvalues <- update_z(hfamily,lb,ub)
##             break
##         }
##     }
##     expect_equal(sum(is.infinite(zvalues)), 0)
## })


## ############################
## ## Test: calculate_scdf() ##
## ############################
## M = 1000
## lb = -5
## ub = 5
## f = dnor
## test_that("The samples created are all with the Boundary)", {
##     library(ars)
##     initialized_sample <- initialize_sample(M, lb, ub, width = 0.5, mod_val = 0)
##     hfamily <- initialized_sample$hfamily
##     zvalues <- initialized_sample$zvalues
##     zcdf <- calculate_scdf(vals = zvalues, hfamily, z = zvalues)
##     expect_equal(zcdf$scdf[1],0)
##     expect_equal(zcdf$scdf[3],1)
## })


## ###########################
## ## Test: create_sample() ##
## ###########################
## M = 1000
## lb =-5
## ub = 5
## f = dnor
## test_that("The samples created are all with the Boundary)", {
##     library(ars)
##     initialized_sample <- initialize_sample(M, lb, ub, width = 0.5, mod_val = 0)
##     hfamily <- initialized_sample$hfamily
##     zvalues <- initialized_sample$zvalues
##     samples <- create_samples(M1 = M^(2/3)*2, hfamily, zvalues)
##     expect_equal(sum(samples<lb), 0)
##     expect_equal(sum(samples>ub), 0)
## })
