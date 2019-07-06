context("setting HmscRandomLevel")

test_that("HmscRandomLevel returns error when no input is given", {
   expect_error(HmscRandomLevel())
})

test_that("HmscRandomLevel returns error when units are specified more than once", {
   expect_error(HmscRandomLevel(units=1:10,sData=data.frame(s1=c(1:10),s2=c(10:1))))
})

test_that("HmscRandomLevel is set properly", {
   expect_equal(HmscRandomLevel(units=1:10)$N,
                10)
   expect_equal(HmscRandomLevel(units=sample(1:5,10,replace=TRUE))$N,
                10)
   expect_equal(levels(HmscRandomLevel(N=10)$pi),
                c('1','2','3','4','5','6','7','8','9','10'))
})

test_that("spatial HmscRandomLevel is set properly", {
   expect_equal(HmscRandomLevel(sData=data.frame(s1=c(1:10),s2=c(10:1)))$sDim,
                2)
   expect_equal(HmscRandomLevel(sData=data.frame(s1=c(1:10),s2=c(10:1)))$N,
                10)
   expect_equal(HmscRandomLevel(sData=data.frame(s1=c(1:10),s2=c(10:1)))$s,
                data.frame(s1=c(1:10),s2=c(10:1)))
   expect_error(HmscRandomLevel(sData=data.frame(s1=c(1:10),s2=c(10:1)), distMat = as.matrix(dist(1:10))))
   expect_equal(HmscRandomLevel(distMat = as.matrix(dist(1:10)))$N,
                10)
   expect_equal(HmscRandomLevel(distMat = as.matrix(dist(1:10)))$distMat,
                as.matrix(dist(1:10)))
})

test_that("covariate dependent HmscRandomLevel is set properly", {
   expect_equal(HmscRandomLevel(xData=data.frame(s1=c(1:10),s2=c(10:1)))$xDim,
                2)
   expect_equal(HmscRandomLevel(xData=data.frame(s1=c(1:10),s2=c(10:1)))$N,
                10)
   expect_equal(HmscRandomLevel(xData=data.frame(s1=c(1:10),s2=c(10:1)))$x,
                data.frame(s1=c(1:10),s2=c(10:1)))
})
