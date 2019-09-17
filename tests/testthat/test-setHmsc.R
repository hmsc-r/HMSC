## ---------------------------------------------------------------------------------------
context("setting species data")
test_that("Hmsc returns error when Y is not a matrix", {
   expect_error(Hmsc(Y=data.frame(s1=1:10),XData=data.frame(x1=1:10)),"Hmsc.setData: Y argument must be a matrix of sampling units times species")
})
## ---------------------------------------------------------------------------------------
context("setting environmental data")
test_that("Hmsc returns error both XData and X are given", {
   expect_error(Hmsc(Y=matrix(1:10),XData=data.frame(x1 = 1:10),X=matrix(1:10)),"Hmsc.setData: only single of XData and X arguments must be specified")
})
test_that("Hmsc returns error when XData is not a data frame or list of data frames", {
   expect_error(Hmsc(Y=matrix(1:10),XData=matrix(1:10)),"Hmsc.setData: XData must either a data.frame or a list of data.frame objects")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),XData=list(s1 = 1:10, s2 = 1:10)),"Hmsc.setData: each element of X list must be a data.frame")
})
test_that("Hmsc returns error when X is not a matrix", {
   expect_error(Hmsc(Y=matrix(1:10),X=data.frame(x1 = 1:10)),"Hmsc.setData: X must either a matrix or a list of matrix objects")
   expect_error(Hmsc(Y=matrix(1:10),X=list(s1 = data.frame(x1 = 1:10))),"Hmsc.setData: each element of X list must be a matrix")
})
test_that("Hmsc returns error when X or XData contains NA values", {
   expect_error(Hmsc(Y=matrix(1:10),X=list(x1 = matrix(c(1:9,NA)))),"all elements of X list must contain no NA values")
   expect_error(Hmsc(Y=matrix(1:10),XData=list(x1 = data.frame(c(1:9,NA)))),"all elements of XData list must contain no NA values")
   expect_error(Hmsc(Y=matrix(1:10),XData=data.frame(c(1:9,NA))),"XData must contain no NA values")
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(c(1:9,NA))),"X must contain no NA values")
})
test_that("Hmsc returns error when the size of X or XData is not correct", {
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=list(s1 = matrix(1:10))),"Hmsc.setData: the length of X list argument must be equal to the number of species")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),XData=list(s1=data.frame(x1 = 1:10))),"the length of XData list argument must be equal to the number of species")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),XData=list(s1=data.frame(x1 = 1:10),s2=data.frame(x1=1:9))),"Hmsc.setData: for each element of XData list the number of rows must be equal to the number of sampling units")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),XData=data.frame(x1 = 1:9)),"Hmsc.setData: the number of rows in XData must be equal to the number of sampling units")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=list(s1=matrix(1:18),s2=matrix(1:10))),"Hmsc.setData: for each element of X list the number of rows must be equal to the number of sampling units")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:9)),"Hmsc.setData: the number of rows in X must be equal to the number of sampling units")
})
test_that("Hmsc returns error when X contains more than one intercept column", {
   Xmat = matrix(rep(1,20),nrow=10,ncol=2)
   colnames(Xmat) = c("Intercept","(Intercept)")
   Xdat = data.frame(Xmat)
   XdatList = list(s1=Xdat,s2=Xdat)
   expect_error(Hmsc(Y=matrix(1:10),X=Xmat),"Hmsc.setData: only one column of X matrix could be named Intercept or (Intercept)",fixed=TRUE)
   expect_error(Hmsc(Y=matrix(1:10),XData=Xdat),"Hmsc.setData: only one column of X matrix could be named Intercept or (Intercept)",fixed=TRUE)
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),XData=XdatList),"Hmsc.setData: only one column of X matrix could be named Intercept or (Intercept)",fixed=TRUE)
})
test_that("Hmsc returns error when intercept column does not only contain 1",{
   Xmat = matrix(rep(1,20),nrow=10,ncol=2)
   Xdat = data.frame(Xmat)
   XdatList = list(s1=Xdat,s2=Xdat)
})
## ---------------------------------------------------------------------------------------
context("setting trait data")
test_that("Hmsc returns error when the wrong combination of Tr, TrData and TrFormula is given", {
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), TrData = data.frame(t1 = 1:2), Tr = matrix(1:2)),"Hmsc.setData: at maximum one of TrData and Tr arguments can be specified")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), TrData = data.frame(t1 = 1:2)),"Hmsc.setData: TrFormula argument must be specified if TrData is provided")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), Tr = data.frame(t1 = 1:2)),"Hmsc.setData: Tr must be a matrix")
})
test_that("Hmsc returns error when trait data does not have the right size",{
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), Tr = matrix(1:3)),"Hmsc.setData: the number of rows in Tr should be equal to number of columns in Y")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), TrData = data.frame(t1=1:3),TrFormula = ~t1),"Hmsc.setData: the number of rows in TrData should be equal to number of columns in Y")
})
test_that("Hmsc returns error when there are unknowns in trait data",{
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), Tr = matrix(c(1,NA))),"Hmsc.setData: Tr parameter must not contain any NA values")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), TrData = data.frame(t1=c(1,NA)),TrFormula = ~t1),"Hmsc.setData: TrData parameter must not contain any NA values")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10), TrData =data.frame(t1 = 1:2), TrFormula = ~ t2),"object 't2' not found")
})
test_that("Hmsc returns error when the intercept columns for Tr is not correctly specified",{
   Traits = matrix(rep(1,4),nrow=2,ncol=2)
   colnames(Traits) = c("Intercept","(Intercept)")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),Tr = Traits),"Hmsc.setData: only one column of Tr matrix could be named Intercept or (Intercept)",fixed=TRUE)
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),TrData = data.frame(Traits),TrFormula = ~ Intercept + (Intercept)),"Hmsc.setData: only one column of Tr matrix could be named Intercept or (Intercept)",fixed=TRUE)
   Traits = matrix(rep(0,2),nrow=2,ncol=1)
   colnames(Traits) = c("Intercept")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),Tr = Traits),"Hmsc.setData: intercept column in Tr matrix must be a column of ones")
})
test_that("Hmsc returns error when intercept is specified in TrFormula",{
   Traits = matrix(rep(1,4),nrow=2,ncol=2)
   colnames(Traits) = c("Intercept","x1")
   expect_error(Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),TrData = data.frame(Traits),TrFormula = ~ Intercept + x1),"Hmsc.setData: only one column of Tr matrix could be named Intercept or (Intercept)",fixed=TRUE)
   Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),TrData = data.frame(Traits),TrFormula = ~ x1)
})
## ---------------------------------------------------------------------------------------
context("setting phylogenetic data")
test_that("Hmsc returns error when both a tree and a correlation matrix are specified",{
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10),phyloTree=TD$phy, C=TD$C),
                "Hmsc.setData: at maximum one of phyloTree and C arguments can be specified")
})
test_that("Hmsc returns error when ",{
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10), C=TD$C),
                "Hmsc.setData: the size of square matrix C must be equal to number of species")
})
## ---------------------------------------------------------------------------------------
context("setting latent structure")
test_that("Hmsc returns error when studydesign does not match random levels",{
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10),ranLevels = TD$rL1),
                "Hmsc.setData: studyDesign is empty, but ranLevels is not")
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10),ranLevels = TD$rL1,studyDesign = data.frame(sample = as.factor(1:10))),
                "Hmsc.setData: studyDesign must contain named columns corresponding to all levels listed in ranLevelsUsed")
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10),ranLevels = list(unit =TD$rL1),studyDesign = data.frame(sample = as.factor(1:10))),
                "Hmsc.setData: studyDesign must contain named columns corresponding to all levels listed in ranLevelsUsed")
   Hmsc(Y=matrix(1:10),X=matrix(1:10),ranLevels = list(sample =TD$rL1),studyDesign = data.frame(sample = as.factor(1:10)))
})
test_that("Hmsc returns error when studydesign does not match Y",{
   expect_error(Hmsc(Y=matrix(1:10),X=matrix(1:10),ranLevels = list(unit =TD$rL1),studyDesign = data.frame(sample = as.factor(1:9))),
                "Hmsc.setData: the number of rows in studyDesign must be equal to number of rows in Y")
})
## ---------------------------------------------------------------------------------------
context("scaling")
test_that("Scaling for Y",{
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10))
   expect_equivalent(colMeans(m$YScaled),c(11/2,31/2))
   expect_equivalent(colMeans(m$Y),c(11/2,31/2))
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),YScale = TRUE)
   expect_equivalent(colMeans(m$YScaled),c(0,0))
   expect_equivalent(colMeans(m$Y),c(11/2,31/2))
   expect_equivalent(round(m$YScalePar),matrix(c(6,3,16,3),ncol=2))
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),YScale = TRUE,distr = c('normal','probit'))
   expect_equivalent(colMeans(m$YScaled),c(0,31/2))
   expect_equivalent(colMeans(m$Y),c(11/2,31/2))
   expect_equivalent(round(m$YScalePar),matrix(c(6,3,0,1),ncol=2))
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),YScale = TRUE,distr = c('poisson','lognormal poisson'))
   expect_equivalent(colMeans(m$YScaled),c(11/2,31/2))
   expect_equivalent(colMeans(m$Y),c(11/2,31/2))
})
test_that("Scaling for X",{
   m = Hmsc(Y=matrix(1:10),X=matrix(1:10))
   expect_equivalent(colMeans(m$X),5.5)
   expect_equivalent(round(colMeans(m$XScaled)),1)
   expect_equivalent(round(m$XScalePar),c(0,7))
   m = Hmsc(Y=matrix(1:10),X=matrix(1:10),XScale = FALSE)
   expect_equivalent(colMeans(m$X),5.5)
   expect_equivalent(colMeans(m$XScaled),5.5)
   m = Hmsc(Y=matrix(1:10),XData=data.frame(x1=1:10),XFormula = ~x1)
   expect_equivalent(colMeans(m$X),c(1,5.5))
   expect_equivalent(colMeans(m$XScaled),c(1,0))
   expect_equivalent(round(m$XScalePar),matrix(c(0,1,6,3),ncol=2))
})
test_that("Scaling for traits",{
   Traits = matrix(c(1,1,2,23),nrow=2,ncol=2)
   colnames(Traits) = c("Intercept","x1")
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),TrData = data.frame(Traits),TrFormula = ~ x1)
   expect_equivalent(m$Tr,Traits)
   expect_equivalent(round(colMeans(m$TrScaled)),c(1,0))
   expect_equivalent(round(m$TrScalePar),matrix(c(0,1,12,15),ncol=2))
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),Tr = Traits)
   expect_equivalent(m$Tr,Traits)
   expect_equivalent(round(colMeans(m$TrScaled)),c(1,0))
   expect_equivalent(round(m$TrScalePar),matrix(c(0,1,12,15),ncol=2))
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:10),Tr = Traits,TrScale = FALSE)
   expect_equivalent(m$TrScaled,Traits)
})
## ---------------------------------------------------------------------------------------
context("Setting data model")
test_that("Family set correctly",{
   m = Hmsc(Y=matrix(1:20,nrow=5,ncol=4),X=matrix(1:5),distr=c('probit','poisson','normal','lognormal poisson'))
   expect_equivalent(m$distr[,1],c(2,3,1,3))
})
test_that("Variance set correctly",{
   m = Hmsc(Y=matrix(1:20,nrow=5,ncol=4),X=matrix(1:5),distr=c('probit','poisson','normal','lognormal poisson'))
   expect_equivalent(m$distr[,2],c(0,0,1,1))
})
test_that("Hmsc returns error when unsuitable data model is given",{
   expect_error(Hmsc(Y=matrix(1:10,nrow=5,ncol=2),X=matrix(1:5),distr=c('probit','logit')),
                "Hmsc.setData: some of the distributions ill defined")
})
## ---------------------------------------------------------------------------------------
context("setting variable selection")

## ---------------------------------------------------------------------------------------
test_that("Hmsc is set correctly",{
   m = Hmsc(Y=TD$Y,
               XData=TD$X,
               XFormula=~x1+x2,
               TrData=TD$Tr,
               TrFormula = ~T1 + T2,
               phyloTree=TD$phy,
               ranLevels=list("sample"=TD$rL2,"plot"=TD$rL1),
               studyDesign = TD$studyDesign,
               distr=c("probit"))
   expect_equal(m$Y,TD$m$Y)
   expect_equal(m$YScaled,TD$m$YScaled)
   expect_equal(m$XData,TD$m$XData)
   expect_equal(m$XScaled,TD$m$XScaled)
   expect_equal(m$C,TD$m$C)
   expect_equal(m$Pi,TD$m$Pi)
   expect_equal(m$studyDesign,TD$m$studyDesign)
   expect_equal(m$TrData,TD$m$TrData)
   expect_equal(m$ranLevels,TD$m$ranLevels)
   expect_equal(m$ranLevelsUsed,TD$m$ranLevelsUsed)
   expect_equal(m$dfPi,TD$m$dfPi)
   expect_equal(m$phyloTree,TD$m$phyloTree)
   expect_equal(m$TrScaled,TD$m$TrScaled)
})
