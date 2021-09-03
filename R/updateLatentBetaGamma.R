#' @importFrom stats rgamma
#'
updateLatentBetaGammaRho = function(GammaLatent, BetaLatent, Varphi, rhoLatent, Tr, C, rL, VC,eC,rhopw){
   ns = nrow(Tr)
   nr = length(rL)
   nt = ncol(Tr)
   TtT = crossprod(Tr)

   for (r in seq_len(nr)) {
      iQLatentList = vector("list", nr)
      if(is.null(C)){
         iQLatentList = rep(list(diag(1,ns)), nf)
      } else{
         iQLatentList = lapply(rhoLatent, function(rho) solve(rho*C+diag(1-rho,ns)))
      }
      varphi = Varphi[[r]]
      betalatent = BetaLatent[[r]]
      gammalatent = GammaLatent[[r]]
      rholatent = rhoLatent[[r]]
      sigma2_gammalatent = rL[[r]]$sigma2_gammalatent
      cns = rL[[r]]$cns
      nf = nrow(betalatent)

      # update zeta by exploiting the varphi decomposition and sampling from TruncNorm
      varphip = varphi
      w0varphi = which(varphi==0)
      probit = pnorm(betalatent[w0varphi])
      varphip[w0varphi]= round( runif(length(w0varphi)) > (1-probit)/(1-probit+probit*(1-cns)) )
      lower = upper = rep(0, nf*ns)
      lower[varphip==0] = -Inf
      upper[varphip==1] = Inf
      zeta = matrix(rtruncnorm(ns*nf, a=lower, b=upper, mean=betalatent, sd=1), nf, ns)

      # update betalatent
      RiU = lapply(vinv, function(x) chol(x+diag(1,ns)) )
      U = lapply(RiU, chol2inv)
      vinvTr = lapply(vinv, function(x) x%*%Tr)
      simnorm = matrix(rnorm(nf*ns), ns,nf)
      betalatent = lapply(seq_len(nf),
                          function(x) {U[[x]]%*%vinvTr[[x]]%*%gammalatent[x,]+zeta[x,] +
                                backsolve(RiU[[x]], simnorm[,x]) } )

      # update gammalatent
      RiU2 = chol(TtT+diag(1/sigma2_gammalatent, nt))
      U2 = chol2inv(RiU2)
      simnorm = matrix(rnorm(nf*nt), nt,nf)
      bs = backsolve(RiU2, simnorm)
      gammalatent  = lapply(seq_len(nf), function(x) {U2%*%Tr%*%betalatent[[x]] + bs[,x]})

      # everything factorizes/parallelises across h=1...nf
      # we use eigendecomposition of C
      VCT = crossprod(VC,Tr)
      ZVC = zeta%*%VC
      eQMat = tcrossprod(matrix(eC),rhopw[,1,drop=FALSE]) + tcrossprod(matrix(1,ns,1),1-rhopw[,1,drop=FALSE])
      # logDetQ  = colSums(log(eQ))
      for(h in seq_len(nf)){
         # we sample (gamma|zeta,rho) by marginalizing out the beta
         eQ = eQMat[,rholatent[h]] #rholatent[h]*eC + (1-rholatent[h])
         iSigma = diag(rep(1/sigma2_gammalatent,nt)) + crossprod((eQ+1)^-0.5 * VCT)
         RiSigma = chol(iSigma)
         tmp = backsolve(RiSigma, crossprod(VCT,((eQ+1)^-1 * ZVC[h,])), transpose=TRUE)
         gammalatent[h,] = backsolve(RiSigma, tmp+rnorm(nt))
         # sample (rho|gamma,zeta).... technically possible to sample (rho|zeta), which would reduce autocorrelation,
         # but it is slightly more numerically complex and tedious/impossible to code loop-free without
         # batched linear algebra operations
         mu = Tr%*%gammalatent[h,]
         eps = zeta[h,] - mu
         VCeps = crossprod(VC,eps)
         eQHatMat = eQMat + 1
         logDetQHat = colSums(log(eQHatMat))
         qF = crossprod(eQHatMat^-1, VCeps^2)
         logLike = log(rhopw[,2]) - 0.5*logDetQHat - 0.5*qF
         logLike = logLike - logSumExp(logLike)
         rholatent[h] = sample(nrow(rhopw),1,prob=exp(logLike))
         # sample (beta|zeta,gamma,rho)
      }

      BetaLatent[[r]] = matrix(unlist(betalatent), nf, ns, byrow = T)
      GammaLatent[[r]] = matrix(unlist(gammalatent), nf, nt, byrow = T)
      rhoLatent[[r]] = rholatent
   }
   return(list(BetaLatent, GammaLatent, rhoLatent))
}
