updateNf = function(eta,lambda,alpha,psi,delta, rL, iter){
   nu = rL$nu
   a1 = rL$a1
   b1 = rL$b1
   a2 = rL$a1
   b2 = rL$b2

   c0 = 1
   c1 = 0.0005
   epsilon = 1e-3                            # threshold limit
   prop = 1.00                   # proportion of redundant elements within columns
   prob = 1/exp(c0 + c1*iter)   # probability of adapting

   if(runif(1) < prob){
      ns = ncol(lambda)
      nf = nrow(lambda)
      np = nrow(eta)

      small = abs(lambda) < epsilon;
      smallProp = apply(small, 1, mean)
      indRedundant = smallProp >= prop
      numRedundant = sum(indRedundant)

      if(nf<rL$nfMax && iter>20 && numRedundant==0 && all(smallProp<0.995)){
         nf = nf+1
         lambdaNew = rbind(lambda, matrix(0,1,ns))
         etaNew = matrix(NA,np,nf)
         etaNew[,1:(nf-1)] = eta
         etaNew[,nf] = rnorm(np)
         alphaNew = c(alpha,1)
         psiNew = matrix(NA,nf,ns)
         psiNew[1:(nf-1),] = psi
         psiNew[nf,] = rgamma(ns,nu/2,nu/2)
         deltaNew = matrix(NA,nf,1)
         deltaNew[1:(nf-1),] = delta
         deltaNew[nf,] = rgamma(1,a2,b2)
         eta = etaNew
         lambda = lambdaNew
         alpha = alphaNew
         psi = psiNew
         delta = deltaNew
      } else if(numRedundant>0 && nf>rL$nfMin){
         indNotRed = setdiff(1:nf, indRedundant)
         nf = length(indNotRed)
         lambda = lambda[indNotRed,,drop=FALSE]
         eta = eta[,indNotRed,drop=FALSE]
         psi = psi[indNotRed,,drop=FALSE]
         delta = delta[indNotRed,,drop=FALSE]
         alpha = alpha[indNotRed]
      }
   }
   return(list(eta=eta, lambda=lambda, alpha=alpha, psi=psi, delta=delta))
}
