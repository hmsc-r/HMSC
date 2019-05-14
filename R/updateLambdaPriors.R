#' @importFrom stats rgamma
#' 
updateLambdaPriors = function(Lambda,Delta, rL){
   nr = length(rL)

   Psi = vector("list", nr)
   for(r in seq_len(nr)){
      nu = rL[[r]]$nu
      a1 = rL[[r]]$a1
      b1 = rL[[r]]$b1
      a2 = rL[[r]]$a2
      b2 = rL[[r]]$b2
      delta = Delta[[r]]
      lambda = Lambda[[r]]
      ns = dim(lambda)[2]
      nf = dim(lambda)[1]
      tau = apply(delta,2,cumprod)
      lambda2 = lambda^2
      aPsi = nu/2 + 0.5

      if(is.matrix(lambda)){
         bPsi = nu/2 + 0.5*lambda2 * matrix(tau,nf,ns)
         psi = matrix(rgamma(nf*ns, aPsi, bPsi), nf, ns)
         M = psi * lambda2
         ad = a1 + 0.5*ns*nf
         bd = b1 + 0.5 * sum(tau * rowSums(M)) / delta[1]
         delta[1] = rgamma(1,ad,bd)
         for(h in 1+seq_len(nf-1)){
            tau = apply(delta,2,cumprod)
            ad = a2 + 0.5*ns*(nf-h+1)
            bd = b2 + 0.5 * sum(tau[-c(1:(h-1)),,drop=FALSE] * apply(M[-c(1:(h-1)),,drop=FALSE],1,sum)) / delta[h]
            delta[h] = rgamma(1,ad,bd)
         }
      } else{
         ncr = dim(lambda)[3]
         bPsi = nu/2 + 0.5*lambda2 * array(tau,c(nf,1,ncr))[,rep(1,ns),,drop=FALSE]
         psi = array(rgamma(nf*ns*ncr, aPsi, bPsi), c(nf,ns,ncr))
         M = psi * lambda2
         ad = a1 + 0.5*ns*nf
         bd = b1 + 0.5 * colSums(tau * apply(M,c(1,3),sum)) / matrix(delta[1,], nf, ncr, byrow=TRUE)
         delta[1,] = rgamma(ncr,ad,bd)
         for(h in 1+seq_len(nf-1)){
            tau = apply(delta,2,cumprod)
            ad = a2 + 0.5*ns*(nf-h+1)
            bd = b2 + 0.5 * colSums(tau[-c(1:(h-1)),,drop=FALSE] * apply(M[-c(1:(h-1)),,,drop=FALSE],c(1,3),sum)) / matrix(delta[h,], nf-h+1, ncr, byrow=TRUE)
            delta[h,] = rgamma(ncr,ad,bd)
         }
      }
      Psi[[r]] = psi
      Delta[[r]] = delta
   }
   return(list(Psi=Psi,Delta=Delta))
}
