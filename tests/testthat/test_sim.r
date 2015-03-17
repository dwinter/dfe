context("Sims work")

w <- rma_normal(n=100, a=0.01, Va=0.001, Ve=0.0001, Ut=5)
                

test_that("we can simulate a normal d.f.e",{
    expect_equal(length(w), 100)
})

