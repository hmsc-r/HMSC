#' @importFrom abind adrop

matmulDiagMat = function(d, V){
   return(matrix(d,length(d),ncol(V)) * V)
}
matmulMatDiag = function(V, d){
   return(V * matrix(d,nrow(V),length(d),byrow=TRUE))
}

recFunBilinearDet = function(node, treeList, X, Y, V, rho, iV, V1, V2, logDetV){
   nChild = treeList[[node]]$n
   parentEdgeLen = treeList[[node]]$parentEdgeLen
   XiSYChild = OneiSXChild = OneiSYChild = OneiSOneChild = vector("list", nChild)
   logDetChild = rep(NA, nChild)
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      if(treeList[[childNode]]$n > 0){
         res = recFunBilinearDet(childNode, treeList, X, Y, V, rho, iV, V1, V2, logDetV)
         XiSYChild[[i]] = res$XiSY
         OneiSXChild[[i]] = res$OneiSX
         OneiSYChild[[i]] = res$OneiSY
         OneiSOneChild[[i]] = res$OneiSOne
         logDetChild[i] = res$logDet
      } else{
         S = treeList[[childNode]]$parentEdgeLen * V1 + V2
         X1 = t(adrop(X[,childNode,,drop=FALSE], 2))
         Y1 = t(adrop(Y[,childNode,,drop=FALSE], 2))
         iSX = solve(S, X1)
         iSY = solve(S, Y1)
         XiSYChild[[i]] = crossprod(X1, iSY)
         OneiSXChild[[i]] = iSX
         OneiSYChild[[i]] = iSY
         RS = chol(S)
         OneiSOneChild[[i]] = chol2inv(RS)
         logDetChild[i] = 2 * sum(log(diag(RS)))
      }
   }
   XiSYSum = Reduce("+", XiSYChild)
   OneiSXSum = Reduce("+", OneiSXChild)
   OneiSYSum = Reduce("+", OneiSYChild)
   OneiSOneSum = Reduce("+", OneiSOneChild)
   logDetSum = sum(logDetChild)
   if(parentEdgeLen == 0)
      return(list(XiSY=XiSYSum, OneiSX=OneiSXSum, OneiSY=OneiSYSum, OneiSOne=OneiSOneSum, logDet=logDetSum))

   W = iV + parentEdgeLen * matmulDiagMat(sqrt(rho), matmulMatDiag(OneiSOneSum, sqrt(rho)))
   RW = chol(W)
   logDet = logDetSum + logDetV + 2*sum(log(diag(RW)))
   iRWT_Drho05_OneiSXSum = backsolve(RW, matmulDiagMat(sqrt(rho),OneiSXSum), transpose=TRUE)
   iRWT_Drho05_OneiSYSum = backsolve(RW, matmulDiagMat(sqrt(rho),OneiSYSum), transpose=TRUE)
   iRWT_Drho05_OneiSOneSum = backsolve(RW, matmulDiagMat(sqrt(rho),OneiSOneSum), transpose=TRUE)
   XiSY = XiSYSum - parentEdgeLen * crossprod(iRWT_Drho05_OneiSXSum, iRWT_Drho05_OneiSYSum)
   OneiSX = OneiSXSum - parentEdgeLen * crossprod(iRWT_Drho05_OneiSOneSum, iRWT_Drho05_OneiSXSum)
   OneiSY = OneiSYSum - parentEdgeLen * crossprod(iRWT_Drho05_OneiSOneSum, iRWT_Drho05_OneiSYSum)
   OneiSOne = OneiSOneSum - parentEdgeLen * crossprod(iRWT_Drho05_OneiSOneSum)
   return(list(XiSY=XiSY, OneiSX=OneiSX, OneiSY=OneiSY, OneiSOne=OneiSOne, logDet=logDet))
}

fastPhyloBilinearDet = function(treeList, X, Y, root, iV, rho){
   RiV = chol(iV)
   V = chol2inv(RiV)
   logDetV = -2*sum(log(diag(RiV)))
   V1 = matmulDiagMat(sqrt(rho), matmulMatDiag(V, sqrt(rho)))
   V2 = matmulDiagMat(sqrt(1-rho), matmulMatDiag(V, sqrt(1-rho)))
   res = recFunBilinearDet(root, treeList, X, Y, V, rho, iV, V1, V2, logDetV)
   return(list(XiSY=res$XiSY, logDet=res$logDet))
}

recFunSampleUp_v1 = function(node, treeList, iV, rho, rho2Mat, XTiDX, XTiDS){
   nc = nrow(iV)
   nChild = treeList[[node]]$n
   edgeLenVec = treeList[[node]]$edgeLen
   iSigmaChild_m = iSigmaHat_m = matrix(NA,nc,nChild)
   iSigmaChild = iSigmaHat = array(NA, c(nc,nc,nChild))
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      if(treeList[[childNode]]$n > 0){
         res = recFunSampleUp(childNode, treeList, iV, rho, rho2Mat, XTiDX, XTiDS)
         iSigmaChild_m[,i] = res$iSigma_m
         iSigmaChild[,,i] = res$iSigma
         treeList = res$treeList
      } else{
         treeList[[childNode]]$iS <- (iSigmaAdded = XTiDX[,,childNode])
         treeList[[childNode]]$iSm <- (iSigmaAdded_beta = XTiDS[,childNode])
         rho2 = rho2Mat[,childNode]
         D2_iSigmaAdded = matmulDiagMat(sqrt(rho2), iSigmaAdded)
         W = iV + matmulMatDiag(D2_iSigmaAdded, sqrt(rho2))
         RW = chol(W)
         iRWT_D2_iSigmaAdded = backsolve(RW, D2_iSigmaAdded, transpose=TRUE)
         iSigmaChild[,,i] = iSigmaAdded - crossprod(iRWT_D2_iSigmaAdded)
         iSigmaChild_m[,i] = iSigmaAdded_beta - crossprod(iRWT_D2_iSigmaAdded, backsolve(RW, sqrt(rho2)*iSigmaAdded_beta, transpose=TRUE))
      }
      D1_iSigmaChild = matmulDiagMat(sqrt(rho), iSigmaChild[,,i])
      W = iV + edgeLenVec[i] * matmulMatDiag(D1_iSigmaChild, sqrt(rho))
      RW = chol(W)
      iRWT_D1_iSigmaChild = backsolve(RW, D1_iSigmaChild, transpose=TRUE)
      iSigmaHat[,,i] =  iSigmaChild[,,i] - edgeLenVec[i] * crossprod(iRWT_D1_iSigmaChild)
      iSigmaHat_m[,i] = iSigmaChild_m[,i] - edgeLenVec[i] * crossprod(iRWT_D1_iSigmaChild, backsolve(RW, sqrt(rho)*iSigmaChild_m[,i], transpose=TRUE))
   }
   iSigmaTilde = rowSums(iSigmaHat, dims=2)
   iSigmaTilde_mTilde = rowSums(iSigmaHat_m)
   treeList[[node]]$iSm <- iSigmaTilde_mTilde
   treeList[[node]]$iS <- iSigmaTilde
   return(list(iSigma_m=iSigmaTilde_mTilde, iSigma=iSigmaTilde, treeList=treeList))
}

recFunSampleUp = function(node, treeList, iV, rho, rho2Mat, XTiDX, XTiDS){
   nc = nrow(iV)
   nChild = treeList[[node]]$n
   edgeLenVec = treeList[[node]]$edgeLen
   iSigmaHat_m = matrix(NA,nc,nChild)
   iSigmaHat = array(NA, c(nc,nc,nChild))
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      if(treeList[[childNode]]$n > 0){
         res = recFunSampleUp(childNode, treeList, iV, rho, rho2Mat, XTiDX, XTiDS)
         iSigmaChild_m = res$iSigma_m
         iSigmaChild = res$iSigma
         treeList = res$treeList
         rhoEff = edgeLenVec[i] * rho
      } else{
         treeList[[childNode]]$iS <- (iSigmaChild = XTiDX[,,childNode])
         treeList[[childNode]]$iSm <- (iSigmaChild_m = XTiDS[,childNode])
         rhoEff = edgeLenVec[i] * rho + rho2Mat[,childNode]
      }
      if(all(diag(iSigmaChild) == Inf)){
         iSigmaHat[,,i] = rhoEff^-1 * iV
         iSigmaHat_m[,i] = iSigmaHat[,,i] %*% iSigmaChild_m
      } else{
         D_iSigmaChild = matmulDiagMat(sqrt(rhoEff), iSigmaChild)
         W = iV + matmulMatDiag(D_iSigmaChild, sqrt(rhoEff))
         RW = chol(W)
         iRWT_D_iSigmaChild = backsolve(RW, D_iSigmaChild, transpose=TRUE)
         iSigmaHat[,,i] = iSigmaChild - crossprod(iRWT_D_iSigmaChild)
         iSigmaHat_m[,i] = iSigmaChild_m - crossprod(iRWT_D_iSigmaChild, backsolve(RW, sqrt(rhoEff)*iSigmaChild_m, transpose=TRUE))
      }
   }
   iSigmaTilde = rowSums(iSigmaHat, dims=2)
   iSigmaTilde_mTilde = rowSums(iSigmaHat_m)
   treeList[[node]]$iSm <- iSigmaTilde_mTilde
   treeList[[node]]$iS <- iSigmaTilde
   return(list(iSigma_m=iSigmaTilde_mTilde, iSigma=iSigmaTilde, treeList=treeList))
}

recFunSampleDown = function(node, treeList, V, rho, rho2Mat, sdMult){
   nc = nrow(V)
   nChild = treeList[[node]]$n
   edgeLenVec = treeList[[node]]$edgeLen
   beta = treeList[[node]]$beta
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      iS = treeList[[childNode]]$iS
      iSm = treeList[[childNode]]$iSm
      V1 = matmulDiagMat(sqrt(rho), matmulMatDiag(V, sqrt(rho)))
      if(treeList[[childNode]]$n > 0){
         W = edgeLenVec[i] * V %*% matmulDiagMat(sqrt(rho), matmulMatDiag(iS, sqrt(rho))) %*% V + V
         RW = chol(W)
         U = edgeLenVec[i] * V1 %*% iS + diag(nc)
         mu = solve(U, beta + edgeLenVec[i] * V1%*%iSm)
         BetaRanPart = edgeLenVec[i]^0.5 * sqrt(rho) * (V %*% backsolve(RW, rnorm(nc)))
         treeList[[childNode]]$beta <- mu + sdMult*BetaRanPart
         treeList = recFunSampleDown(childNode, treeList, V, rho, rho2Mat, sdMult)
      } else{
         rho2 = rho2Mat[,childNode]
         W = edgeLenVec[i] * V1 + matmulDiagMat(sqrt(rho2), matmulMatDiag(V, sqrt(rho2)))
         RW = chol(W)
         iSigma = iS + chol2inv(RW)
         iSigma_mu = iSm + backsolve(RW, backsolve(RW, beta, transpose=TRUE))
         RiSigma = chol(iSigma)
         BetaRanPart = backsolve(RiSigma, rnorm(nc))
         treeList[[childNode]]$beta <- backsolve(RiSigma, backsolve(RiSigma, iSigma_mu, transpose=TRUE)) + sdMult*BetaRanPart
      }
   }
   return(treeList)
}

recFunSample = function(treeList, root, V, iV, rho, rho2Mat, XTiDX, XTiDS, sdMult=1){
   ns = ncol(XTiDS)
   nc = nrow(V)
   treeListTemp = treeList
   for(i in 1:length(treeListTemp)){
      treeListTemp[[i]]$iS = NULL
      treeListTemp[[i]]$iSm = NULL
      treeListTemp[[i]]$beta = NULL
   }
   resUp = recFunSampleUp(root, treeListTemp, iV, rho, rho2Mat, XTiDX, XTiDS)
   treeListTemp = resUp$treeList
   treeListTemp[[root]]$beta = matrix(0, nc, 1)
   treeListTemp = recFunSampleDown(root, treeListTemp, V, rho, rho2Mat, sdMult)
   Beta = matrix(NA, nc, ns)
   for(j in 1:ns){
      Beta[,j] = treeListTemp[[j]]$beta
   }
   return(Beta)
}




recFunScalarUp = function(node, treeList, rho, XTiDX, XTiDS){
   nc = nrow(XTiDS)
   nChild = treeList[[node]]$n
   edgeLenVec = treeList[[node]]$edgeLen
   iSigmaHat_m = matrix(NA, nc, nChild)
   iSigmaHat = rep(NA, nChild)
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      if(treeList[[childNode]]$n > 0){
         res = recFunScalarUp(childNode, treeList, rho, XTiDX, XTiDS)
         iSigmaChild_m = res$iSigma_m
         iSigmaChild = res$iSigma
         treeList = res$treeList
         rhoEff = edgeLenVec[i] * rho
      } else{
         treeList[[childNode]]$iS <- (iSigmaChild = XTiDX[childNode])
         treeList[[childNode]]$iSm <- (iSigmaChild_m = XTiDS[,childNode])
         rhoEff = edgeLenVec[i] * rho + (1 - rho)
      }
      if(iSigmaChild == Inf){
         iSigmaHat[i] = rhoEff^-1
         iSigmaHat_m[,i] = iSigmaHat[i] * iSigmaChild_m
      } else if(iSigmaChild == 0){
         iSigmaHat[i] = 0
         iSigmaHat_m[,i] = iSigmaChild_m
      } else{
         iSigmaHat[i] = (iSigmaChild^-1 + rhoEff)^-1
         iSigmaHat_m[,i] = (iSigmaChild_m / iSigmaChild) * iSigmaHat[i]
      }
   }
   iSigmaTilde = sum(iSigmaHat)
   iSigmaTilde_mTilde = rowSums(iSigmaHat_m)
   treeList[[node]]$iSm <- iSigmaTilde_mTilde
   treeList[[node]]$iS <- iSigmaTilde
   return(list(iSigma_m=iSigmaTilde_mTilde, iSigma=iSigmaTilde, treeList=treeList))
}

recFunScalarDown = function(node, treeList, rho, sdMult){
   # nc = nrow(V)
   nChild = treeList[[node]]$n
   edgeLenVec = treeList[[node]]$edgeLen
   beta = treeList[[node]]$beta
   for(i in seq_len(nChild)){
      childNode = treeList[[node]]$child[i]
      iS = treeList[[childNode]]$iS
      iSm = treeList[[childNode]]$iSm
      if(treeList[[childNode]]$n > 0){
         rhoEff = edgeLenVec[i] * rho
      } else{
         rhoEff = edgeLenVec[i] * rho + (1 - rho)
      }
      if(iS == Inf){
         mu = iSm
      } else{
         isigma2 = rhoEff^-1 + iS
         mu = (iSm + beta / rhoEff) / isigma2
      }
      
      treeList[[childNode]]$beta <- mu #+ sdMult*BetaRanPart
      if(treeList[[childNode]]$n > 0){
         treeList = recFunScalarDown(childNode, treeList, rho, sdMult)
      }
   }
   return(treeList)
}

recFunScalar = function(treeList, root, rho, XTiDX, XTiDS, sdMult=1){
   ns = ncol(XTiDS)
   nc = nrow(XTiDS)
   treeListTemp = treeList
   for(i in 1:length(treeListTemp)){
      treeListTemp[[i]]$iS = NULL
      treeListTemp[[i]]$iSm = NULL
      treeListTemp[[i]]$beta = NULL
   }
   resUp = recFunScalarUp(root, treeListTemp, rho, XTiDX, XTiDS)
   treeListTemp = resUp$treeList
   treeListTemp[[root]]$beta = matrix(0, nc, 1)
   treeListTemp = recFunScalarDown(root, treeListTemp, rho, sdMult)
   Beta = matrix(NA, nc, ns)
   for(j in 1:ns){
      Beta[,j] = treeListTemp[[j]]$beta
   }
   return(Beta)
}



recFunScalarUp_opt = function(node, envTree, rho, XTiDX, XTiDS){
   nodeEnv = envTree[[node]]
   nChild = nodeEnv$n
   edgeLenVec = nodeEnv$edgeLen
   
   # OPTIMIZATION: Instead of allocating iSigmaHat_m (matrix) and iSigmaHat (vector) 
   # just to sum them later, we initialize running accumulators. This prevents 
   # thousands of small matrix allocations during recursion.
   iSigmaTilde = 0
   iSigmaTilde_mTilde = numeric(nrow(XTiDS))
   
   for(i in seq_len(nChild)){
      childIdx = nodeEnv$child[i]
      childEnv = envTree[[childIdx]]
      
      if(childEnv$n > 0){
         res = recFunScalarUp_opt(childIdx, envTree, rho, XTiDX, XTiDS)
         iSigmaChild_m = res$iSigma_m
         iSigmaChild = res$iSigma
         rhoEff = edgeLenVec[i] * rho
      } else {
         # Leaf node: modify environment in-place
         iSigmaChild = XTiDX[childIdx]
         iSigmaChild_m = XTiDS[, childIdx]
         
         childEnv$iS <- iSigmaChild
         childEnv$iSm <- iSigmaChild_m
         
         rhoEff = edgeLenVec[i] * rho + (1 - rho)
      }
      
      # Calculate current child's contribution
      if(iSigmaChild == Inf){
         iSigmaHat_i = rhoEff^-1
         iSigmaHat_m_i = iSigmaHat_i * iSigmaChild_m
      } else if(iSigmaChild == 0){
         iSigmaHat_i = 0
         iSigmaHat_m_i = iSigmaChild_m
      } else {
         iSigmaHat_i = (iSigmaChild^-1 + rhoEff)^-1
         iSigmaHat_m_i = (iSigmaChild_m / iSigmaChild) * iSigmaHat_i
      }
      
      # Accumulate running sums directly
      iSigmaTilde = iSigmaTilde + iSigmaHat_i
      iSigmaTilde_mTilde = iSigmaTilde_mTilde + iSigmaHat_m_i
   }
   
   # Update parent environment in-place
   nodeEnv$iSm <- iSigmaTilde_mTilde
   nodeEnv$iS <- iSigmaTilde
   
   return(list(iSigma_m = iSigmaTilde_mTilde, iSigma = iSigmaTilde))
}

recFunScalarDown_opt = function(node, envTree, rho, sdMult){
   nodeEnv = envTree[[node]]
   nChild = nodeEnv$n
   edgeLenVec = nodeEnv$edgeLen
   beta = nodeEnv$beta
   
   for(i in seq_len(nChild)){
      childIdx = nodeEnv$child[i]
      childEnv = envTree[[childIdx]]
      
      iS = childEnv$iS
      iSm = childEnv$iSm
      
      if(childEnv$n > 0){
         rhoEff = edgeLenVec[i] * rho
      } else {
         rhoEff = edgeLenVec[i] * rho + (1 - rho)
      }
      
      if(iS == Inf){
         mu = iSm
      } else {
         isigma2 = rhoEff^-1 + iS
         mu = (iSm + beta / rhoEff) / isigma2
      }
      
      # Update beta in-place
      childEnv$beta <- mu 
      
      if(childEnv$n > 0){
         recFunScalarDown_opt(childIdx, envTree, rho, sdMult)
      }
   }
   return(invisible(NULL))
}

recFunScalar_opt = function(treeList, root, rho, XTiDX, XTiDS, sdMult=1){
   ns = ncol(XTiDS)
   nc = nrow(XTiDS)
   
   # Convert to environments (pass-by-reference).
   # OPTIMIZATION: Skipped the loop that NULLs out variables. It is unnecessary 
   # because the Up/Down passes will strictly overwrite them anyway.
   envTree = lapply(treeList, function(x) list2env(x, parent = emptyenv()))
   
   # Upward Pass
   recFunScalarUp_opt(root, envTree, rho, XTiDX, XTiDS)
   
   # Downward Pass Setup
   envTree[[root]]$beta = matrix(0, nc, 1)
   
   # Downward Pass
   recFunScalarDown_opt(root, envTree, rho, sdMult)
   
   # Extract final results
   # OPTIMIZATION: Initializing with NA_real_ avoids implicit type coercion 
   # from logical NA to double when populating the matrix.
   Beta = matrix(NA_real_, nc, ns)
   for(j in seq_len(ns)){
      Beta[,j] = envTree[[j]]$beta
   }
   
   return(Beta)
}