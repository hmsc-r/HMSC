#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom tensorflow tf
#'
updateEta = function(Z,Beta,iD,Eta,Lambda,Alpha, rLPar, X,Pi,dfPi,rL, tfCompFlag=FALSE){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   bandMatrix = function(A){
      d1 = dim(A)[1]
      d2 = dim(A)[2]
      ind1 = rep(1:(d1*d2), each=d2)
      ind2 = rep(rep((0:(d2-1))*d1,d1) + rep(1:d1,each=d2), d2)
      bMat = sparseMatrix(ind1, ind2, x=as.vector(aperm(A,c(3,1,2))))
   }

   switch(class(X)[1L],
          matrix = {
             LFix = X%*%Beta
          },
          list = {
             LFix = matrix(NA,ny,ns)
             for(j in 1:ns)
                LFix[,j] = X[[j]]%*%Beta[,j]
          }
   )
   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,ny,ns)
         for(k in 1:rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
      }
   }
   for(r in seq_len(nr)){
      rnames=rownames(Eta[[r]])
      if(nr > 1){
         S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      } else{
         S = Z - LFix
      }
      iDS = iD*S
      iDS[is.na(Z)] = 0
      lambda = Lambda[[r]]
      nf = dim(lambda)[1]
      lPi = Pi[,r]
      ldfPi = dfPi[,r]
      randEps = matrix(rnorm(np[r]*nf),np[r],nf)

      if(tfCompFlag==TRUE){
         tfla = tf$linalg
         ic = function(...) as.integer(c(...))
         if(rL[[r]]$xDim == 0){
            LamiDLam = tf$scatter_nd(matrix(ic(lPi)), tfla$einsum("hj,ij,gj->ihg",lambda,iD,lambda), tf$constant(ic(np[r],nf,nf)))
            PiDSLambda = tf$matmul(tf$scatter_nd(matrix(ic(lPi)), iDS, tf$constant(ic(np[r],ns))), lambda, transpose_b=TRUE)
         } else{
            lambdaLocal = tf$einsum("ik,hjk->ihj", rL[[r]]$x[unLdfPi,], lambda)
            LamiDLam = tf$scatter_nd(matrix(ic(lPi)), tfla$einsum("ihj,ij,igj->ihg",lambdaLocal,iD,lambdaLocal),
                                     tf$constant(ic(np[r],nf,nf)))
            PiDSLambda = tf$scatter_nd(matrix(ic(lPi)), tf$matmul(iDS, lambdaLocal, transpose_b=TRUE), tf$constant(ic(np[r],ns)))
         }
      } else{
         LamiDLam = array(0,c(np[r],nf,nf))
         PiDSLambda = matrix(0,np[r],nf)
         for(q in 1:np[r]){
            rows = which(lPi==q)
            if(rL[[r]]$xDim == 0){
               for(p in seq_len(length(rows)))
                  LamiDLam[q,,] = LamiDLam[q,,] + tcrossprod(lambda*matrix(sqrt(iD[rows[p],]),nf,ns,byrow=TRUE))
               PiDSLambda[q,] = colSums(tcrossprod(iDS[rows,,drop=FALSE], lambda))
            } else{
               ncr = rL[[r]]$xDim
               lambdaLocal = rowSums(lambda * array(unlist(rep(rL[[r]]$x[unLdfPi[q],],each=nf*ns)), c(nf,ns,ncr)), dims=2)
               for(p in seq_len(length(rows)))
                  LamiDLam[q,,] = LamiDLam[q,,] + tcrossprod(lambdaLocal*matrix(sqrt(iD[rows[p],]),nf,ns,byrow=TRUE))
               PiDSLambda[q,] = colSums(tcrossprod(iDS[rows,,drop=FALSE], lambdaLocal))
            }
         }
      }

      if(rL[[r]]$sDim == 0){
         eta = matrix(NA,np[r],nf)
         if(tfCompFlag==TRUE){
            iV = diag(nf) + LamiDLam
            LiV = tfla$cholesky(iV)
            m = tfla$cholesky_solve(LiV, tf$expand_dims(PiDSLambda, ic(-1)))
            res = m + tfla$triangular_solve(LiV, tf$expand_dims(randEps,ic(-1)), adjoint=TRUE)
            eta = tf$squeeze(res, ic(-1))$numpy()
         } else{
            for(q in 1:np[r]){
               rows = which(lPi==q)
               iV = diag(nf) + LamiDLam[q,,]
               RiV = chol(iV)
               V = chol2inv(RiV)
               mu = PiDSLambda[q,] %*% V
               eta[q,] = mu + backsolve(RiV,randEps[q,])
            }
         }
      } else{
         eta = matrix(0,np[r],nf)
         alpha = Alpha[[r]]
         iWg = rLPar[[r]]$iWg
         if(rL[[r]]$spatialMethod == "Full"){
            iWs = bdiag(lapply(seq_len(nf), function(x) iWg[,,alpha[x]]))
            LamiDLamBandMat = bandMatrix(LamiDLam)
            iUEta = iWs + LamiDLamBandMat
            R = chol(iUEta)
            tmp2 = backsolve(R, as.vector(PiDSLambda), transpose=TRUE) + as.vector(randEps)
            feta = backsolve(R, tmp2);
            eta = matrix(feta,np[r],nf);
         } else if(rL[[r]]$spatialMethod == "NNGP"){
            iWs = sparseMatrix(c(),c(),dims=c(np[r]*nf,np[r]*nf))
            for(h in seq_len(nf))
               iWs = iWs + kronecker(iWg[[alpha[h]]], Diagonal(x=c(rep(0,h-1),1,rep(0,nf-h))))
            LamiDLamBlockMat = bdiag(lapply(split(LamiDLam, 1:np[r]), function(x) matrix(x,nf,nf)))
            iUEta = iWs + LamiDLamBlockMat
            R = chol(iUEta)
            tmp2 = backsolve(R, as.vector(t(PiDSLambda)), transpose=TRUE) + as.vector(t(randEps))
            feta = backsolve(R, tmp2)
            eta = matrix(feta,np[r],nf,byrow=TRUE)
         } else if(rL[[r]]$spatialMethod == "GPP"){
            # here iD and idD correspond to observation noise and GPP approximation parts respectively
            idDg = rLPar[[r]]$idDg
            idDW12g = rLPar[[r]]$idDW12g
            Fg = rLPar[[r]]$Fg
            nK = nrow(Fg)
            idD = idDg[,alpha]
            Fmat = matrix(0,nrow=(nK*nf),ncol=(nK*nf))
            idD1W12 = matrix(0,nrow=(np[r]*nf),ncol=(nK*nf))
            for(h in 1:nf){
               Fmat[(h-1)*nK+(1:nK), (h-1)*nK+(1:nK)] = Fg[,,alpha[h]]
               idD1W12[(h-1)*np[r]+(1:np[r]), (h-1)*nK+(1:nK)] = idDW12g[,,alpha[h]]
            }

            iAst = LiAst = array(NA,c(np[r],nf,nf))
            for(i in 1:np[r]){
               Ael = cholLamiDLam[i,,]+diag(idD[i,])
               RAel = chol(Ael)
               iAst[i,,] = chol2inv(RAel)
               LiAst[i,,] = t(chol(iAst[i,,]))
            }
            iA = bandMatrix(B)
            LiA = bandMatrix(LB)
            iAidD1W12 = iA %*% idD1W12
            H = Fmat - t(idD1W12)%*%iAidD1W12
            RH = chol(as.matrix(H))
            iRH = solve(RH)

            mu1 = iA%*%as.vector(PiDSLambda)
            tmp1 = iAidD1W12 %*% iRH
            mu2 = tmp1%*%(Matrix::t(tmp1)%*%as.vector(PiDSLambda))

            etaR = LiA%*%rnorm(np[r]*nf) + tmp1%*%rnorm(nK*nf)
            eta = matrix(mu1+mu2+etaR,ncol=nf,nrow=np[r])
         }
      }
      rownames(eta)=rnames
      Eta[[r]] = eta
      if(r < nr){
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,ny,ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
         }
      }
   }
   return(Eta)
}

