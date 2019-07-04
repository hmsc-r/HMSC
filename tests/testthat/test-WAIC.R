context("compute WAIC")

test_that("HmscRandomLevel is set properly", {
   expect_equal(round(computeWAIC(TD$m),1),
                1.8)
})
