#' @importFrom stats rgamma rbinom
#'
#'
# NEW version that updates only continuous prior Psi, Delta
# Psi is the local continuous scale
# Delta is the column continuous scale
# Vartheta is the remaining part of the column scale that can be
# either mixture of normal (CUSP) or product of delta (MGP, in this case Vartheta = NULL)

updateLambdaContPriors = function(LambdaTilde, Delta, Vartheta, rL){
   nr = length(rL)
   Psi = vector("list", nr)
   for (r in seq_len(nr)) {
      nu = rL[[r]]$nu
      a1 = rL[[r]]$a1
      b1 = rL[[r]]$b1
      a2 = rL[[r]]$a2
      b2 = rL[[r]]$b2
      delta_1 = Delta[[r]]   # LS: changed the name in delta_1, consistently with our notation
      lambdaTilde = LambdaTilde[[r]]
      vartheta = Vartheta[[r]]  # LS: added vartheta
      ns = dim(lambdaTilde)[2]
      nf = dim(lambdaTilde)[1]

      if(rL[[r]]$progShrinkType=="MGP"){
         tau_1 = apply(delta_1, 2, cumprod)     # classic MGP
      } else if(rL[[r]]$progShrinkType=="CUSP"){
         tau_1 = delta_1 / vartheta               # LS: Tau^(-1) is (delta X vartheta)^(-1)
      }
      lambdaTilde2 = lambdaTilde^2
      aPsi = nu/2 + 0.5
      if (is.matrix(lambdaTilde)) {
         bPsi = nu/2 + 0.5 * lambdaTilde2 * matrix(tau_1, nf, ns)
         psi_1 = matrix(rgamma(nf * ns, aPsi, bPsi), nf, ns)  # LS: changed the name in psi_1, consistently with our notation
         M = psi_1 * lambdaTilde2
         rowSumM = rowSums(M)
         aVec = c(a1,rep(a2,nf-1))
         bVec = c(b1,rep(b2,nf-1))
         if(rL[[r]]$progShrinkType=="MGP"){
            for(h in seq_len(nf)) {
               tau_1 = apply(delta_1, 2, cumprod)
               ad = aVec[h] + 0.5 * ns*(nf-h+1)
               bd = bVec[h] + 0.5 * sum(tau_1[h:nf, ,drop=FALSE] * rowSumM[h:nf]) / delta_1[h]
               delta_1[h] = rgamma(1, ad, bd)
            }
         } else if(rL[[r]]$progShrinkType=="CUSP"){
            adVec = aVec + 0.5 * ns
            bdVec = aVec + 0.5 * rowSumM / vartheta  # LS: added division by vartheta
            delta_1 = rgamma(nf, adVec, bdVec)
         }
      }
      else {
         #GT: lets omit this more complicated case for now
         stop("Capacity for covariate-dependent associations is currently disabled")
         ncr = dim(lambdaTilde)[3]
         bPsi = nu/2 + 0.5 * lambda2 * array(tau, c(nf,1,ncr))[,rep(1,ns),,drop = FALSE]
         psi = array(rgamma(nf*ns*ncr, aPsi, bPsi), c(nf,ns,ncr))
         M = psi * lambdaTilde2

         if(rL[[r]]$progShrinkType=="MGP"){
            ad = a1 + 0.5 * ns * nf
            bd = b1 + 0.5 * colSums(tau * apply(M, c(1, 3), sum))/matrix(delta[1,], nf, ncr, byrow=TRUE)
         } else if(rL[[r]]$progShrinkType=="CUSP"){
            ad = a1 + 0.5 * ns
            bd = b1 + 0.5 * apply(M, c(1, 3), sum)[1,]
         }
         delta[1, ] = rgamma(ncr, ad, bd)

         for (h in 1+seq_len(nf-1)) {
            if(rL[[r]]$progShrinkType=="MGP"){
               tau = apply(delta, 2, cumprod)
               ad = a2 + 0.5 * ns * (nf - h + 1)
               bd = b2 + 0.5 * colSums(tau[h:nf,,drop = FALSE] * apply(M[-c(1:(h - 1)), , ,
                                                                         drop = FALSE], c(1, 3), sum))/matrix(delta[h,
                                                                         ], nf - h + 1, ncr, byrow = TRUE)
            } else if(rL[[r]]$progShrinkType=="CUSP"){
               ad = a2 + 0.5 * ns
               bd = b2 + 0.5 * apply(M, c(1, 3), sum)[h,]
            }
            delta[h,] = rgamma(ncr, ad, bd)
         }
      }
      Psi[[r]] = psi_1
      Delta[[r]] = delta_1
   }
   return(list(Psi=Psi, Delta=Delta))
}
