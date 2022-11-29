# @title updateBetaLambda
#
# @description updates beta lambda
#
#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix
#'
updateBetaLambda = function(Z,Gamma,iV,iD,Eta,Psi,Delta, iQ, X,Tr,Pi,dfPi,C,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   S = Z
   Lambda = vector("list", nr)
   if(nr > 0){
      EtaFull = vector("list", nr)
      nf = rep(NA,nr)
      ncr = rep(NA,nr)
      for(r in seq_len(nr)){
         if(rL[[r]]$xDim == 0){
            EtaFull[[r]] = Eta[[r]][Pi[,r],]
         } else{
            EtaFull[[r]] = vector("list", rL[[r]]$xDim)
            for(k in 1:rL[[r]]$xDim)
               EtaFull[[r]][[k]] = Eta[[r]][Pi[,r],] * rL[[r]]$x[as.character(dfPi[,r]),k]
         }
         nf[r] = ncol(Eta[[r]])
         ncr[r] = max(rL[[r]]$xDim, 1)
      }
      nfSum = sum(nf*ncr)
      EtaSt = matrix(unlist(EtaFull),ny,nfSum)
      if(class(X)[1L] == "matrix"){
         XEta = cbind(X, EtaSt)
      } else if(class(X)[1L] == "list"){
         XEta = lapply(X, cbind, EtaSt)
      }
      PsiT = vector("list",nr)
      for(r in 1:nr){
         if(rL[[r]]$xDim == 0){
            PsiT[[r]] = t(Psi[[r]])
         } else{
            PsiT[[r]] = aperm(Psi[[r]], c(2,1,3))
         }
      }
      psiSt = matrix(unlist(PsiT),nfSum,ns,byrow=TRUE)
      Tau = lapply(Delta, function(a) apply(a,2,cumprod))
      tauSt = matrix(unlist(Tau),nfSum,1)
      priorLambda = psiSt*matrix(tauSt,nfSum,ns)
   } else{
      nf = c()
      ncr = c()
      nfSum = 0
      XEta = X
      priorLambda = matrix(numeric(0),0,ns)
   }

   Mu = rbind(tcrossprod(Gamma,Tr), matrix(0,nfSum,ns))
   iDS = iD*S
   iDS[is.na(Z)] = 0
   XEtaTiDXEtaList = vector("list", ns)
   if(class(X)[1L] == "matrix"){
      XTiDS = crossprod(XEta,iDS)
      for(j in 1:ns){
         XEtaTiDXEtaList[[j]] = crossprod(sqrt(iD[,j])*XEta)
      }
   } else if(class(X)[1L] == "list"){
      XTiDS = matrix(NA,nc+nfSum,ns)
      for(j in 1:ns){
         XTiDS[,j] = crossprod(XEta[[j]],iDS[,j])
         XEtaTiDXEtaList[[j]] = crossprod(sqrt(iD[,j])*XEta[[j]])
      }
   }

   if(is.null(C)){ # no phylogeny information
      BetaLambda = matrix(NA, nc+nfSum, ns)
      randEps = matrix(rnorm((nc+nfSum)*ns), nc+nfSum, ns)
      for(j in 1:ns){
         P = bdiag(iV, diag(priorLambda[,j]))
         iU = P + XEtaTiDXEtaList[[j]]
         RiU = chol(iU)
         m1 = backsolve(RiU, P%*%Mu[,j] + XTiDS[,j], transpose=TRUE)
         BetaLambda[,j] = backsolve(RiU, m1 + randEps[,j])
      }
   } else{ # available phylogeny information
      P = bdiag(kronecker(iV,iQ), Diagonal(x=t(priorLambda)))
      tmpMat = Reduce(rbind, XEtaTiDXEtaList)
      ind1 = rep(rep(1:ns,each=nc+nfSum)+ns*rep(0:(nc+nfSum-1),ns), nc+nfSum)
      ind2 = rep(1:((nc+nfSum)*ns), each=nc+nfSum)
      XEtaTiDXEtaMat = sparseMatrix(ind1, ind2, x=as.vector(tmpMat))
      RiU = chol(as.matrix(XEtaTiDXEtaMat) + P)
      m1 = backsolve(RiU, P%*%as.vector(t(Mu)) + as.vector(t(XTiDS)), transpose=TRUE)
      BetaLambda = matrix(backsolve(RiU, m1 + rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
   }
   Beta = BetaLambda[1:nc,,drop=FALSE]
   nfCumSum = c(0,cumsum(nf*ncr)) + nc
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         Lambda[[r]] = BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE]
      } else
         Lambda[[r]] = aperm(array(BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE], c(nf[r],ncr[r],ns)),c(1,3,2))
   }
   return(list(Beta=Beta, Lambda=Lambda))
}

