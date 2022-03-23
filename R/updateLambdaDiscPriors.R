#' @importFrom stats rgamma rmultinom
#' @importFrom ggdist dstudent_t

updateLambdaDiscPriors = function(Z, Beta,iSigma, Eta,Lambda,LambdaTilde,Delta, BetaLatent,Varphi,Vartheta,TildeU,W, X,Pi,rL){
   # require(mc2d)    # vectorize density of multinomial likelihood
   # require(matrixStats)
   nr = length(rL)
   ny = nrow(Z)
   sigma05 = 1/sqrt(iSigma)

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

   Lambda = vector("list", nr)
   for(r in seq_len(nr)){
      rnames=rownames(Eta[[r]])
      if(nr > 1){
         S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      } else{
         S = Z - LFix
      }
      epsilon = S

      xi = rL[[r]]$xi
      cns = rL[[r]]$cns
      # LS: added vartheta_inf and delta_h hyperparameters
      a1 = rL[[r]]$a1
      b1 = rL[[r]]$b1
      a2 = rL[[r]]$a2
      b2 = rL[[r]]$b2
      vartheta_inf = rL[[r]]$vartheta_inf
      eta = Eta[[r]]
      lambdaTilde = LambdaTilde[[r]]
      delta_1 = Delta[[r]]  #LS: added Delta, it is used to compute tau.
      betaLatent = BetaLatent[[r]]
      varphi = Varphi[[r]]
      vartheta = Vartheta[[r]]
      tildeU = TildeU[[r]]    # LS: added tildeU among function parameters
      w = W[[r]]
      ns = dim(lambdaTilde)[2]
      nf = dim(lambdaTilde)[1]

      probit = pnorm(betaLatent)*cns
      logPrior1 = log(probit)
      logPrior0 = log(1-probit)

      if(rL[[r]]$xDim==0){
         # if(rL[[r]]$progShrinkType=="CUSP"){
         #    actFact = which(vartheta==1)
         #    if(nf > length(actFact)){
         #       inactFact = setdiff(seq_len(nf), actFact)
         #       # sumlog = matrix( apply( cbind(as.vector(logprobit[inactFact,]),
         #       #                               as.vector(logprobit_1[inactFact,])), 1, logSumExp), inactFact,ns)
         #       # varphi[inactFact, ] = round( matrix(runif(ns*length(inactFact)),inactFact,ns) <
         #       #                                   exp(logprobit_1[inactFact,]-sumlog) )
         #       varphi[inactFact,] = rbinom(length(inactFact)*ns, 1, probit[inactFact,])
         #    }
         # } else if(rL[[r]]$progShrinkType=="MGP"){
         #    actFact = seq_len(nf)
         # }
         # LS: commented above if statement.
         # LS: we can now consider all active factors also in CUSP to sample local scale
         actFact = seq_len(nf)

         #GT_old: maybe we can sample for several h from joint distribution. The complexity grows as (ns*nf/|h| * 2^|h|), so for e.g. 10 randomly selected factors shall be OKish
         # However, coding this efficiently seems rather nontrivial, so leave the idea for later.
         if(rL[[r]]$structuredShrinkage==TRUE){
            for(h in actFact){
               # update eta %*% lambda
               etalambda = eta[,actFact]%*%(lambdaTilde[actFact,]*varphi[actFact,])
               # NormalMean0 = etalambda - t(varphi[h,]*lambdaTilde[h, ]*matrix(eta[,h], ns, ny, byrow = T))
               NormalMean0 = etalambda - eta[,h,drop=FALSE]%*%(lambdaTilde[h,,drop=FALSE]*varphi[h,,drop=FALSE])
               logLikeNorm0 = colSums(dnorm(epsilon, mean=NormalMean0, sd=matrix(sigma05,ny,ns,byrow=TRUE), log=TRUE))
               NormalMean1 = NormalMean0 + eta[,h,drop=FALSE]%*%lambdaTilde[h,,drop=FALSE]
               logLikeNorm1 = colSums(dnorm(epsilon, mean=NormalMean1, sd=matrix(sigma05,ny,ns,byrow=TRUE), log=TRUE))

               lp0 = logPrior0[h,] + logLikeNorm0
               lp1 = logPrior1[h,] + logLikeNorm1
               sumlog = apply(cbind(lp0, lp1),1, logSumExp)
               varphi[h, ] = rbinom(ns, 1, exp(lp1-sumlog))
            }
         }else{
            varphi = matrix(1,nf,ns)
         }

         # we udate vartheta if it is not null
         if(rL[[r]]$progShrinkType=="CUSP"){
            u = rep(NA,nf)
            # LS: now lambdaTilde includes vartheta
            vartheta_sqrt = sqrt(vartheta)
            lambdaTilde_scaled =  lambdaTilde / vartheta_sqrt
            for(h in seq_len(nf)){
               # LS: now lambdaTilde includes vartheta
               etalambda = eta%*%(lambdaTilde*varphi)
               tmp_scaled_factor_h =  eta[,h,drop=FALSE]%*%(lambdaTilde_scaled[h,,drop=FALSE]*varphi[h,,drop=FALSE])

               # LS: sampling the square root of vartheta
               NormalVariance = (sum(tmp_scaled_factor_h^2 * matrix(iSigma,ny,ns,byrow=TRUE)) + 1/tildeU[h])^(-1)
               mean_adjust_h = (epsilon - etalambda + eta[,h,drop=FALSE]%*%(lambdaTilde[h,,drop=FALSE]*varphi[h,,drop=FALSE])) * matrix(iSigma,ny,ns,byrow=TRUE) #GT: added division by sigma
               NormalMean = NormalVariance * sum(tmp_scaled_factor_h*mean_adjust_h)
               # print(c(NormalMean, NormalVariance))
               vartheta_sqrt[h] = rnorm(1, NormalMean, sqrt(NormalVariance)) #GT: added sqrt, since rnorm operates with sd
               # LS: update lambda tilde
               lambdaTilde = lambdaTilde_scaled * vartheta_sqrt
            }
            # LS: update the augmented data u and tildeU
            vartheta = vartheta_sqrt^2
            tau_sqrt = vartheta_sqrt / sqrt(delta_1)
            # LS: faster way to compute probabilities
            # loglik0_1 = dstudent_t(tau_sqrt[1], df=2*a_1, mu=0, sigma=sqrt(vartheta_inf*b_1/a_1), log=T)
            # loglik1_1 = dstudent_t(tau_sqrt[1], df=2*a_1, mu=0, sigma=sqrt(b_1/a_1), log=T)
            # loglik0 = c(loglik0_1, dstudent_t(tau_sqrt[-1],df=2*a_2, mu=0, sigma=sqrt(vartheta_inf*b_2/a_2), log=T))
            # loglik1 = c(loglik1_1, dstudent_t(tau_sqrt[-1], df=2*a_2, mu=0, sigma=sqrt(b_2/a_2), log=T))
            loglik0 = dstudent_t(tau_sqrt, df=2*c(a1,rep(a2,nf-1)), mu=0, sigma=sqrt(vartheta_inf*c(b1/a1,rep(b2/a2,nf-1))), log=T)
            loglik1 = dstudent_t(tau_sqrt, df=2*c(a1,rep(a2,nf-1)), mu=0, sigma=sqrt(c(b1/a1,rep(b2/a2,nf-1))), log=T)
            loglik_mat = matrix(rep(loglik0,nf), nf, nf, byrow=T)
            loglik_uppmat = matrix(rep(loglik1,nf), nf, nf, byrow=T)
            loglik_mat[upper.tri(loglik_mat)] = loglik_uppmat[upper.tri(loglik_uppmat)]
            logProb = loglik_mat + matrix(rep(log(w),nf), nf, nf, byrow=T)
            for(h in seq_len(nf)){
               prob1 = exp(logProb[h,] - logSumExp(logProb[h,]))
               u[h] = which(rmultinom(n=1, size=1, prob=prob1)==1)
               tildeU[h] = vartheta_inf + as.numeric(u[h]>h)*(1-vartheta_inf)
            }

            tmp = table(c(u,1:nf)) - 1
            v = rbeta(nf, shape1=1+tmp, shape2=xi+nf-cumsum(tmp))
            v[nf] = 1
            w = v*c(1,cumprod(1-v[-nf]))
         }
      }else{
         stop("Capacity for covariate-dependent associations is currently disabled")
      }

      if(rL[[r]]$progShrinkType=="CUSP"){
         Vartheta[[r]] = vartheta
         TildeU[[r]] = tildeU
         W[[r]] = w
         LambdaTilde[[r]] = lambdaTilde  # LS: saved also lambdaTilde
      }
      Varphi[[r]] = varphi
      Lambda[[r]] = LambdaTilde[[r]]*Varphi[[r]]
   }
   return(list(Lambda=Lambda, Varphi=Varphi, Vartheta=Vartheta, W=W, TildeU=TildeU))
}
