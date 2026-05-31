context("Sanity recovery for all updaters")

# Skip these computationally expensive tests unless explicitly requested
skip_if_not(Sys.getenv("RUN_HMSC_SANITY_TESTS") == "true",
            message = "Skipping expensive sanity recovery tests. Set RUN_HMSC_SANITY_TESTS='true' to run.")

# Helper function to compute index recovery success
test_index_recovery = function(samples, true_val) {
  res = table(samples)
  as.numeric(names(res)[which.max(res)])
}

# 1. updateRho
test_that("updateRho recovers rhoInd correctly (General and Boundary)", {
  set.seed(42)
  ns = 15
  nc = 2
  rhopw = cbind(seq(0, 1, length.out=11), rep(1/11, 11))
  phy = ape::rtree(ns)
  phy$tip.label = paste0("sp", 1:ns)
  phy = ape::compute.brlen(phy)
  C = ape::vcv(phy, corr=TRUE)

  true_rho_ind_b = 10 
  rho_b = rhopw[true_rho_ind_b, 1]
  Gamma = matrix(0, nc, 1)
  Tr = matrix(1, ns, 1)
  Beta_b = matrix(NA, nc, ns)
  for(k in 1:nc){
     Beta_b[k, ] = MASS::mvrnorm(1, mu=rep(0, ns), Sigma=rho_b*C + (1-rho_b)*diag(ns))
  }
  Y_dummy = matrix(0, 10, ns)
  colnames(Y_dummy) = phy$tip.label
  m = Hmsc(Y=Y_dummy, X=matrix(rnorm(10*nc), 10, nc), phyloTree=phy)
  m$rhopw = rhopw
  dp = computeDataParameters(m)
  rho_samples_b = numeric(50)
  for(i in 1:50){
    rho_samples_b[i] = updateRho(Beta=Beta_b, Gamma=Gamma, iV=diag(nc),
                                 RQg=dp$RQg, detQg=dp$detQg, Tr=Tr, rhopw=rhopw)
  }
  expect_true(abs(test_index_recovery(rho_samples_b, true_rho_ind_b) - true_rho_ind_b) <= 2)
})

# 2. updateAlpha
test_that("updateAlpha recovers AlphaInd correctly (General and Boundary)", {
  set.seed(42)
  np = 80
  nf = 1
  coords = matrix(runif(np*2), np, 2)
  rownames(coords) = 1:np
  rL = HmscRandomLevel(sData = coords)
  rL = setPriors(rL)
  true_alpha_ind_b = 7
  alpha_b = rL$alphapw[true_alpha_ind_b, 1]
  Sigma_b = exp(-as.matrix(dist(coords))/alpha_b)
  eta_b = matrix(MASS::mvrnorm(1, mu=rep(0, np), Sigma=Sigma_b), np, nf)
  m = Hmsc(Y=matrix(0, np, 2), X=matrix(1, np, 1), ranLevels=list("unit"=rL), studyDesign=data.frame(unit=as.factor(1:np)))
  dp = computeDataParameters(m)
  alpha_samples_b = numeric(50)
  for(i in 1:50){
    res = updateAlpha(Eta=list(eta_b), rL=m$rL, rLPar=dp$rLPar)
    alpha_samples_b[i] = res[[1]][1]
  }
  expect_true(abs(test_index_recovery(alpha_samples_b, true_alpha_ind_b) - true_alpha_ind_b) <= 3)
})

# 3. updateEta
test_that("updateEta recovers latent factors", {
  set.seed(42)
  ny = 40
  ns = 5
  nf = 1
  true_eta = matrix(rnorm(ny), ny, nf)
  Lambda = list(matrix(rnorm(nf*ns), nf, ns))
  Z = true_eta %*% Lambda[[1]] + matrix(rnorm(ny*ns, sd=0.05), ny, ns)
  rL = HmscRandomLevel(N=ny)
  m = Hmsc(Y=matrix(0, ny, ns), X=matrix(1, ny, 1), ranLevels=list("u"=rL), studyDesign=data.frame(u=as.factor(1:ny)))
  dp = computeDataParameters(m)
  eta_res = updateEta(Y=m$Y, Z=Z, Beta=matrix(0, 1, ns), iSigma=rep(400, ns), Eta=list(matrix(0, ny, nf)),
                      Lambda=Lambda, AlphaInd=list(c(1)), rLPar=dp$rLPar, Loff=matrix(0, ny, ns),
                      X=m$X, Pi=m$Pi, dfPi=m$dfPi, rL=m$rL)
  expect_true(cor(as.vector(eta_res[[1]]), as.vector(true_eta)) > 0.8)
})

# 4. updateBetaLambda
test_that("updateBetaLambda recovers Beta correctly", {
  set.seed(42)
  ny = 50
  ns = 4
  nc = 2
  X = matrix(rnorm(ny*nc), ny, nc)
  true_beta = matrix(rnorm(nc*ns), nc, ns)
  Z = X %*% true_beta + matrix(rnorm(ny*ns, sd=0.05), ny, ns)
  beta_lambda = updateBetaLambda(Y=matrix(0, ny, ns), Z=Z, Gamma=matrix(0, nc, 1), iV=diag(nc),
                                 iSigma=rep(400, ns), Eta=list(), Psi=list(), Delta=list(),
                                 rhoInd=1, iQ=diag(ns), Loff=matrix(0, ny, ns), X=X, Tr=matrix(1, ns, 1),
                                 Pi=matrix(nrow=ny, ncol=0), dfPi=data.frame(row.names=1:ny), C=diag(ns), rL=list())
  expect_true(cor(as.vector(beta_lambda$Beta), as.vector(true_beta)) > 0.9)
})

# 5. updateInvSigma
test_that("updateInvSigma recovers residual variance", {
  set.seed(42)
  ny = 200
  ns = 5
  true_sigma = 0.4
  E = matrix(rnorm(ny*ns, sd=sqrt(true_sigma)), ny, ns)
  iSigma = updateInvSigma(Y=matrix(0, ny, ns), Z=E, Beta=matrix(0, 1, ns), iSigma=rep(1, ns),
                          Eta=list(), Lambda=list(), distr=matrix(c(1, 1), ns, 2, byrow=TRUE),
                          Loff=matrix(0, ny, ns), X=matrix(0, ny, 1), Pi=matrix(nrow=ny, ncol=0),
                          dfPi=data.frame(row.names=1:ny), rL=list(), aSigma=rep(0.1, ns), bSigma=rep(0.1, ns))
  expect_true(abs(mean(1/iSigma) - true_sigma) < 0.1)
})

# 6. updateZ
test_that("updateZ respects truncation for probit", {
  set.seed(42)
  ny = 60
  ns = 2
  Y = matrix(sample(0:1, ny*ns, replace=TRUE), ny, ns)
  Mu = matrix(rnorm(ny*ns), ny, ns)
  Z = updateZ(Y=Y, Z=matrix(0, ny, ns), Beta=matrix(0, 1, ns), iSigma=rep(1, ns),
              Eta=list(), Lambda=list(), Loff=Mu, X=matrix(0, ny, 1), Pi=matrix(nrow=ny, ncol=0),
              dfPi=data.frame(row.names=1:ny), distr=matrix(c(2, 0), ns, 2, byrow=TRUE), rL=list(), ind=NULL)
  expect_true(all(Z[Y == 1] > 0))
  expect_true(all(Z[Y == 0] < 0))
})

# 7. updateGammaV
test_that("updateGammaV recovers Gamma and V", {
  set.seed(42)
  ns = 80
  nc = 2
  nt = 2
  Tr = matrix(rnorm(ns*nt), ns, nt)
  true_gamma = matrix(rnorm(nc*nt), nc, nt)
  true_V = diag(nc) * 0.15
  Beta = true_gamma %*% t(Tr) + t(MASS::mvrnorm(ns, mu=rep(0, nc), Sigma=true_V))
  res = updateGammaV(Beta=Beta, Gamma=matrix(0, nc, nt), iV=diag(nc), rhoInd=1,
                     Tr=Tr, C=diag(ns), iQg=array(diag(ns), c(ns, ns, 1)),
                     RQg=array(diag(ns), c(ns, ns, 1)), mGamma=matrix(0, nc*nt, 1),
                     iUGamma=diag(nc*nt), V0=diag(nc), f0=nc+1)
  expect_true(cor(as.vector(res$Gamma), as.vector(true_gamma)) > 0.9)
})

# 8. updateLambdaPriors (Psi and Delta)
test_that("updateLambdaPriors recovery", {
  set.seed(42)
  nf = 5
  ns = 30
  Lambda = list(matrix(rnorm(nf*ns, sd=0.05), nf, ns))
  Lambda[[1]][1, ] = rnorm(ns, sd=2)
  rL = HmscRandomLevel(N=30)
  rL = setPriors(rL, nfMax=nf, nfMin=1)
  res = updateLambdaPriors(Lambda=Lambda, Delta=list(matrix(1, nf, 1)), rL=list(rL))
  expect_true(res$Delta[[1]][1] < mean(res$Delta[[1]][2:nf]))
})

# 9. updateNf
test_that("updateNf handles redundancy", {
  set.seed(42)
  nf = 10
  ns = 5
  np = 5
  lambda = matrix(0, nf, ns)
  eta = matrix(rnorm(np*nf), np, nf)
  alphaInd = rep(1, nf)
  psi = matrix(1, nf, ns)
  delta = matrix(100, nf, 1)
  rL = HmscRandomLevel(N=np)
  rL$nfMin = 1
  rL$nfMax = 20
  rL$nu = 1
  for(i in 1:20){
     res = updateNf(eta=eta, lambda=lambda, alphaInd=alphaInd, psi=psi, delta=delta, rL=rL, iter=i)
     if(ncol(res$eta) < nf) break
  }
  expect_true(ncol(res$eta) < nf)
})

# 10. updateBetaSel
test_that("updateBetaSel handles variable selection", {
  set.seed(42)
  ns = 4
  ny = 30
  Z = matrix(rnorm(ny*ns), ny, ns)
  X = list()
  for(j in 1:ns) X[[j]] = matrix(rnorm(ny*2), ny, 2)
  XSelect = list(list(covGroup=2, spGroup=rep(1, ns), q=rep(0.5, 1)))
  res = updateBetaSel(Z=Z, XSelect=XSelect, BetaSel=list(rep(TRUE, 1)), Beta=matrix(0, 2, ns),
                      iSigma=rep(1, ns), Lambda=list(), Eta=list(), Loff=matrix(0, ny, ns),
                      X1=X, Pi=matrix(nrow=ny, ncol=0), dfPi=data.frame(row.names=1:ny), rL=list())
  expect_true(is.list(res$BetaSel))
  expect_true(is.logical(res$BetaSel[[1]]))
})

# 11. updateGamma2
test_that("updateGamma2 recovers Gamma correctly", {
  set.seed(42)
  ny = 100
  ns = 10
  nc = 2
  nt = 1
  X = matrix(rnorm(ny*nc), ny, nc)
  Tr = matrix(1, ns, nt)
  true_gamma = matrix(c(1, -1), nc, nt)
  Beta = true_gamma %*% t(Tr) + matrix(rnorm(nc*ns, sd=0.1), nc, ns)
  Z = X %*% Beta + matrix(rnorm(ny*ns, sd=0.1), ny, ns)
  
  res = updateGamma2(Z=Z, iV=diag(nc), iSigma=rep(1, ns), Eta=list(), Lambda=list(), 
                     Loff=NULL, X=X, Pi=matrix(nrow=ny, ncol=0), dfPi=data.frame(row.names=1:ny),
                     Tr=Tr, C=NULL, rL=list(), iQg=NULL, mGamma=matrix(0, nc*nt, 1), iUGamma=diag(nc*nt))
  
  expect_true(cor(as.vector(res), as.vector(true_gamma)) > 0.95)
})

# 12. updateGammaEta
test_that("updateGammaEta execution and basic recovery", {
  set.seed(42)
  ny = 50
  ns = 5
  nc = 1
  nt = 1
  nf = 1
  X = matrix(1, ny, nc)
  Tr = matrix(1, ns, nt)
  true_gamma = matrix(0.5, nc, nt)
  true_eta = matrix(rnorm(ny), ny, nf)
  Lambda = list(matrix(1, nf, ns))
  Z = X %*% true_gamma %*% t(Tr) + true_eta %*% Lambda[[1]] + matrix(rnorm(ny*ns, sd=0.1), ny, ns)
  
  rL = HmscRandomLevel(N=ny)
  dp = computeDataParameters(Hmsc(Y=matrix(0, ny, ns), X=X))
  
  res = updateGammaEta(Z=Z, Gamma=matrix(0, nc, nt), V=diag(nc), iV=diag(nc), id=rep(1, ns),
                       Eta=list(matrix(0, ny, nf)), Lambda=Lambda, AlphaInd=list(1),
                       Loff=NULL, X=X, Tr=Tr, Pi=matrix(1:ny, ny, 1), dfPi=data.frame(u=as.factor(1:ny)),
                       rL=list(rL), rLPar=dp$rLPar, Q=diag(ns), iQ=diag(ns), RQ=diag(ns),
                       mGamma=matrix(0, nc*nt, 1), U=diag(nc*nt), iU=diag(nc*nt))
  
  expect_true(is.matrix(res$Gamma))
  expect_true(cor(as.vector(res$Eta[[1]]), as.vector(true_eta)) > 0.7)
})

# 13. updateLatentLoadingOrder
test_that("updateLatentLoadingOrder executes correctly", {
  set.seed(42)
  nf = 3
  ns = 10
  np = 10
  eta = matrix(rnorm(np*nf), np, nf)
  lambda = matrix(rnorm(nf*ns), nf, ns)
  alphaInd = c(1, 2, 3)
  delta = matrix(c(1, 1.5, 2), nf, 1)
  rL = HmscRandomLevel(N=np)
  rL = setPriors(rL)
  
  res = updateLatentLoadingOrder(eta=eta, lambda=lambda, alphaInd=alphaInd, delta=delta, rL=rL)
  expect_equal(dim(res$eta), c(np, nf))
  expect_equal(dim(res$lambda), c(nf, ns))
  expect_equal(length(res$alpha), nf)
})

# 14. updatewRRR and updatewRRRPriors
test_that("updatewRRR and updatewRRRPriors execution", {
  set.seed(42)
  ny = 40
  ns = 4
  ncNRRR = 1
  ncRRR = 2
  ncORRR = 2
  X1A = matrix(1, ny, ncNRRR)
  XRRR = matrix(rnorm(ny*ncORRR), ny, ncORRR)
  Beta = matrix(rnorm((ncNRRR+ncRRR)*ns), ncNRRR+ncRRR, ns)
  Z = matrix(rnorm(ny*ns), ny, ns)
  
  PsiRRR = matrix(1, ncRRR, ncORRR)
  DeltaRRR = matrix(1, ncRRR, 1)
  
  rL = list()
  res = updatewRRR(Z=Z, Beta=Beta, iSigma=rep(1, ns), Eta=list(), Lambda=list(),
                   Loff=NULL, X1A=X1A, XRRR=XRRR, Pi=matrix(nrow=ny, ncol=0),
                   dfPi=data.frame(row.names=1:ny), rL=rL, PsiRRR=PsiRRR, DeltaRRR=DeltaRRR)
  
  expect_equal(dim(res$wRRR), c(ncRRR, ncORRR))
  
  res_p = updatewRRRPriors(wRRR=res$wRRR, Delta=DeltaRRR, nu=1, a1=1, b1=1, a2=1, b2=1)
  expect_equal(dim(res_p$Psi), c(ncRRR, ncORRR))
})
