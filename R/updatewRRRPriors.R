#' @importFrom stats rgamma
#' 
updatewRRRPriors = function(wRRR,Delta,nu,a1,b1,a2,b2){
   delta = Delta
   lambda = wRRR
   ns = dim(lambda)[2]
   nf = dim(lambda)[1]
   tau = apply(delta,2,cumprod)
   lambda2 = lambda^2
   aPsi = nu/2 + 0.5

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
   Psi = psi
   Delta = delta
   return(list(Psi=Psi,Delta=Delta))
}
