#' @importFrom stats rgamma
#'
updateLatentBetaGammaRho = function(GammaLatent,BetaLatent,Varphi,rhoLatent, Tr,phyloPar,rL, rhopw){
   ns = nrow(Tr)
   nr = length(rL)
   nt = ncol(Tr)
   TtT = crossprod(Tr)
   eC = phyloPar$eC
   VC = phyloPar$VC

   for (r in seq_len(nr)) {
      # iQLatentList = vector("list", nr)
      # if(is.null(C)){
      #    iQLatentList = rep(list(diag(1,ns)), nf)
      # } else{
      #    iQLatentList = lapply(rhoLatent, function(rho) solve(rho*C+diag(1-rho,ns)))
      # }
      varphi = Varphi[[r]]
      betaLatent = BetaLatent[[r]]
      gammaLatent = GammaLatent[[r]]
      rholatent = rhoLatent[[r]]
      sigma2_gammaLatent = rL[[r]]$sigma2_gammaLatent
      cns = rL[[r]]$cns
      nf = nrow(betaLatent)

      # update zeta by exploiting the varphi decomposition and sampling from TruncNorm
      varphip = varphi
      w0varphi = which(varphi==0)
      probit = pnorm(betaLatent[w0varphi])
      varphip[w0varphi] = round( runif(length(w0varphi)) > (1-probit)/(1-probit+probit*(1-cns)) )
      lower = upper = rep(0, nf*ns)
      lower[varphip==0] = -Inf
      upper[varphip==1] = Inf
      zeta = matrix(rtruncnorm(nf*ns, a=lower, b=upper, mean=betaLatent, sd=1), nf, ns)

      # # update betaLatent
      # RiU = lapply(vinv, function(x) chol(x+diag(1,ns)) )
      # U = lapply(RiU, chol2inv)
      # vinvTr = lapply(vinv, function(x) x%*%Tr)
      # simnorm = matrix(rnorm(nf*ns), ns,nf)
      # betaLatent = lapply(seq_len(nf),
      #                     function(x) {U[[x]]%*%vinvTr[[x]]%*%gammaLatent[x,]+zeta[x,] +
      #                           backsolve(RiU[[x]], simnorm[,x]) } )
      #
      # # update gammaLatent
      # RiU2 = chol(TtT+diag(1/sigma2_gammaLatent, nt))
      # U2 = chol2inv(RiU2)
      # simnorm = matrix(rnorm(nf*nt), nt,nf)
      # bs = backsolve(RiU2, simnorm)
      # gammaLatent  = lapply(seq_len(nf), function(x) {U2%*%Tr%*%betaLatent[[x]] + bs[,x]})

      # everything factorizes/parallelises across h=1...nf
      # we use eigendecomposition of C
      VCT = crossprod(VC,Tr)
      ZVC = zeta%*%VC
      eQMat = tcrossprod(matrix(eC),rhopw[,1,drop=FALSE]) + tcrossprod(matrix(1,ns,1),1-rhopw[,1,drop=FALSE])
      eQHatMat = eQMat + 1
      logDetQHat = colSums(log(eQHatMat))
      for(h in seq_len(nf)){
         # we sample (gamma|zeta,rho) by marginalizing out the beta
         eQ = eQMat[,rholatent[h]]
         iSigma = diag(rep(sigma2_gammaLatent^-1,nt)) + crossprod((eQ+1)^-0.5 * VCT)
         RiSigma = chol(iSigma)
         tmp = backsolve(RiSigma, crossprod(VCT,((eQ+1)^-1 * ZVC[h,])), transpose=TRUE)
         gammaLatent[h,] = backsolve(RiSigma, tmp+rnorm(nt))

         # sample (rho|gamma,zeta).... technically possible to sample (rho|zeta), which would reduce autocorrelation,
         # but it is slightly more numerically complex and tedious/impossible to code loop-free without
         # batched linear algebra operations
         mu = Tr%*%gammaLatent[h,]
         eps = zeta[h,] - mu
         VCeps = crossprod(VC,eps)
         qF = crossprod(eQHatMat^-1, VCeps^2)
         logLike = log(rhopw[,2]) - 0.5*logDetQHat - 0.5*qF
         logLike = logLike - logSumExp(logLike)
         rholatent[h] = sample(nrow(rhopw),1,prob=exp(logLike))

         # sample (beta|zeta,gamma,rho)
         eQ = eQMat[,rholatent[h]] #rholatent[h]*eC + (1-rholatent[h])
         eK = (eQ^-1 + 1^-1)^-1
         v = eQ^-1*(VCT%*%gammaLatent[h,]) + ZVC[h,]
         betaLatent[h,] = VC%*%(eK*v + rnorm(ns,0,sqrt(eK)))
      }

      BetaLatent[[r]] = betaLatent
      GammaLatent[[r]] = gammaLatent
      rhoLatent[[r]] = rholatent
   }
   return(list(BetaLatent=BetaLatent, GammaLatent=GammaLatent, rhoLatent=rhoLatent))
}
