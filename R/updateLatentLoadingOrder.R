#' @importFrom stats runif rnorm dgamma dt
#' @importFrom abind abind
#' @importFrom pracma meshgrid
#' @importFrom matrixStats logSumExp
updateLatentLoadingOrder = function(eta,lambda,alpha,delta, rL){
   #!!!!!! Current implementation will not work for covariate-dependnet latent loadings!!!!!!!!!!
   if(rL$xDim <= 1){
      nu = rL$nu
      a1 = rL$a1
      b1 = rL$b1
      a2 = rL$a2
      b2 = rL$b2
      ns = dim(lambda)[2]
      nf = dim(lambda)[1]

      # print(dim(eta))
      # print(dim(lambda))
      # print(dim(alpha))
      # print(dim(delta))

      # these values shall better depend on the a/b parameters
      deltaLogMin = -1
      deltaLogMax = 5
      N = 100
      aArray = array(rep(c(a1,rep(a2,nf-1)), each=N^2), c(N,N,nf))
      bArray = array(rep(c(b1,rep(b2,nf-1)), each=N^2), c(N,N,nf))

      deltaLogGrid = abind(rev(meshgrid(seq(deltaLogMin,deltaLogMax,length.out=N))), along=3)
      tauLogGrid = aperm(apply(deltaLogGrid,c(1,2),cumsum), c(2,3,1))
      deltaLog = log(delta[,1])
      for(h in 1:(nf-1)){
         tauLog = cumsum(deltaLog)
         tauLogMat = array(rep(tauLog,each=N^2),c(N,N,nf))
         tauLogMat[,,h:(h+1)] = tauLog[h]-deltaLog[h]+tauLogGrid
         deltaLogMat = abind(tauLogMat[,,1], aperm(array(apply(tauLogMat,c(1,2),diff),c(nf-1,N,N)), c(2,3,1)), along=3)
         tmp = array(exp(0.5*tauLogMat[,,h:(h+1)]), c(N,N,2,ns))
         Lambda0Mat = array(rep(lambda[h:(h+1),],each=N^2), c(N,N,2,ns)) * tmp
         Lambda0Mat_swap = array(rep(lambda[(h+1):h,],each=N^2), c(N,N,2,ns)) * tmp
         deltaMat = exp(deltaLogMat)
         log_prior_delta = rowSums(dgamma(deltaMat,aArray,bArray,log=TRUE), dims=2) + rowSums(deltaLogMat, dims=2)
         log_prior_Lambda = rowSums(dt(Lambda0Mat,nu,log=TRUE), dims=2)
         log_prior_Lambda_swap = rowSums(dt(Lambda0Mat_swap,nu,log=TRUE), dims=2)
         logSpaceWarp = 0.5*ns*deltaLogGrid[,,1] + 0.5*ns*(deltaLogGrid[,,1]+deltaLogGrid[,,2])
         log_prior_loading = log_prior_delta + log_prior_Lambda + logSpaceWarp
         log_prior_loading_swap = log_prior_delta + log_prior_Lambda_swap + logSpaceWarp
         log_prob = logSumExp(log_prior_loading) - 2*log(N)
         log_prob_swap = logSumExp(log_prior_loading_swap) - 2*log(N)
         swapRatio = 1 / (exp(log_prob-log_prob_swap) + 1)
         if(runif(1) <= swapRatio){
            cat(sprintf("swapping %d %d, prob=%.4f\n", h, h+1, swapRatio))
            ind = sample(0:(N^2-1), 1, prob=exp(log_prior_loading_swap-max(log_prior_loading_swap)))
            ind1 = ind%%N + 1
            ind2 = ind%/%N + 1
            lambda[h:(h+1),] = lambda[(h+1):h,]
            eta[,h:(h+1)] = eta[,(h+1):h]
            alpha[h:(h+1)] = alpha[(h+1):h]
         } else{
            ind = sample(0:(N^2-1), 1, prob=exp(log_prior_loading-max(log_prior_loading)))
            ind1 = ind%%N + 1
            ind2 = ind%/%N + 1
         }
         deltaLog = deltaLogMat[ind1,ind2,]
      }
      delta = matrix(exp(deltaLog), nf,1)
   }
   return(list(eta=eta, lambda=lambda, alpha=alpha, delta=delta))
}
