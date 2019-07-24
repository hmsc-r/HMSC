context("compute WAIC")

test_that("WAIC gives correct output", {
   expect_equal(round(computeWAIC(TD$m),1),
                0.8)
})
