## ---------------------------------------------------------------------------------------
context("setting priors for random levels")
test_that("Hmsc returns error when nfMin > nfMax", {
   expect_error(setPriors(TD$rL1, nfMin=10, nfMax=5),
                "HmscRandomLevel.setPriors: nfMin must be not greater than nfMax")
})
test_that("Hmsc returns error when alpha prior is incorrect", {
   expect_error(setPriors(TD$rL2, alphapw = matrix(c(0,1,1,2),nrow=2,ncol=2)),
                "HmscRandomLevel.setPriors: prior for spatial scale was given, but not spatial coordinates were specified")
   expect_error(setPriors(TD$rL1, alphapw = matrix(c(0,1,1,2))),
                "HmscRandomLevel.setPriors: alphapw must be a matrix with two columns")
})
test_that("Priors are correct", {
   expect_equal(setPriors(TD$rL1,nfMax=2, nfMin=2),TD$rL1)
   expect_equal(setPriors(TD$rL2,nfMax=2, nfMin=2),TD$rL2)
})
## ---------------------------------------------------------------------------------------
context("setting priors for fixed effects")
test_that("Hmsc returns error when phylo prior is incorrect", {
   m = Hmsc(Y = matrix(1:10),X = matrix(1:10))
   expect_error(setPriors(m, rhopw = matrix(c(0,1,1,2),nrow=2,ncol=2)),
                "HMSC.setPriors: prior for phylogeny given, but no phylogenic relationship matrix was specified")
   expect_error(setPriors(TD$m, rhopw = matrix(c(0,1,1,2))),
                "HMSC.setPriors: rhopw must be a matrix with two columns")
})
test_that("Priors are correct", {
   expect_equal(setPriors(TD$m),TD$m)
})

