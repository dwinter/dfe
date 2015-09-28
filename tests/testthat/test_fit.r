context("Building fitting functions")

test_that("We can generate box constrainsts", {
    starts  <- list(shape=1, rate=10)
    both <- list(shape=1e-8, rate=1e-8)
    one <- list(shape=1e-8)
    none <- list()
    expect_equal(process_constraints(starts, both), c(1e-8, 1e-8))
    expect_equal(process_constraints(starts, one),  c(1e-8, -Inf))
    expect_null(process_constraints(starts, none, FALSE))
})

test_that( "We can make custom fitting fxns", {
   f <- make_dfe_fitting_fxn(dma_gamma, c("shape", "rate", "Ve", "Ut"))
   expect_is(f, "function")
   expect_equal(names(formals(f)), c("obs", "start", "fixed", "verbose"))
})
