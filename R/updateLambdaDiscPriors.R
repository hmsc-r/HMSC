#' @importFrom stats rgamma

updateLambdaDiscPriors = function(Z, Beta,iSigma, Eta,Lambda,LambdaTilde,BetaLatent,Varphi,Vartheta,W, X,Pi,rL){
   #GT: I think multinomial is covered in "stats" and we do not need multinormal density as covariance matrices are diagonal.
   # require(mc2d)    # vectorize density of multinomial likelihood
   # require(matrixStats)
   nr = length(rL)
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
      eta = Eta[[r]]
      lambdaTilde = LambdaTilde[[r]]
      betaLatent = BetaLatent[[r]]
      varphi = Varphi[[r]]
      vartheta = Vartheta[[r]]
      w = W[[r]]
      ns = dim(lambdaTilde)[2]
      nf = dim(lambdaTilde)[1]

      probit = pnorm(betaLatent)*cns
      logPrior1 = log(probit)
      logPrior0 = log(1-probit)

      if(rL[[r]]$xDim==0){
         if(rL[[r]]$progShrinkType=="CUSP"){
            actFact = which(vartheta==1)
            if(nf > length(actFact)){
               inactFact = setdiff(seq_len(nf), actFact)
               # sumlog = matrix( apply( cbind(as.vector(logprobit[inactFact,]),
               #                               as.vector(logprobit_1[inactFact,])), 1, logSumExp), inactFact,ns)
               # varphi[inactFact, ] = round( matrix(runif(ns*length(inactFact)),inactFact,ns) <
               #                                   exp(logprobit_1[inactFact,]-sumlog) )
               varphi[inactFact, ] = rbinom(length(inactFact)*ns, 1, probit[inactFact,])
            }
         } else if(rL[[r]]$progShrinkType=="MGP"){
            actFact = seq_len(nf)
         }
         #GT: maybe we can sample for several h from joint distribution. The complexity grows as (ns*nf/|h| * 2^|h|), so for e.g. 10 randomly selected factors shall be OKish
         # However, coding this efficiently seems rather nontrivial, so leave the idea for later.
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

         # we udate vartheta if it is not null
         if(rL[[r]]$progShrinkType=="CUSP"){
            u = rep(NA,nf)
            for(h in seq_len(nf)){
               etalambda = eta%*%(lambdaTilde*varphi*vartheta)
               NormalMean0 = etalambda - eta[,h,drop=FALSE]%*%(lambdaTilde[h,,drop=FALSE]*varphi[h,,drop=FALSE]*vartheta[h])
               logLikeNorm0 = sum(dnorm(epsilon, mean=NormalMean0, sd=matrix(sigma05,ny,ns,byrow=TRUE), log=TRUE))
               NormalMean1 = NormalMean0 + eta[,h,drop=FALSE]%*%(lambdaTilde[h,,drop=FALSE]*varphi[h,,drop=FALSE])
               logLikeNorm1 = sum(dnorm(epsilon, mean=NormalMean1, sd=matrix(sigma05,ny,ns,byrow=TRUE), log=TRUE))

               logProb1 = c(rep(logLikeNorm0,h), rep(logLikeNorm1,nf-h)) + log(w)
               prob1 = exp(logProb1 - logSumExp(logProb1))
               u[h] = which(rmultinom(n=1, size=1, prob=prob1)==1)
               vartheta[h] = as.numeric(u[h]>h)
            }
            tmp = table(c(u,1:nf)) - 1
            v = rbeta(nf, shape1=1+tmp, shape2=xi+nf-cumsum(tmp))
            v[nf] = 1
            w = v*c(1,cumprod(1-v[-nf])) #GT: could it more numerically stable to operate with log(w), log(v), log(1-v)?
         }
      }else{
         stop("Capacity for covariate-dependent associations is currently disabled")
      }

      Varphi[[r]] = varphi
      Lambda[[r]] = LambdaTilde[[r]]*Varphi[[r]] #GT: we also need to update the Lambda here to be used on other conditional updaters
      if(rL[[r]]$progShrinkType=="CUSP"){
         Vartheta[[r]] = vartheta
         W[[r]] = w
         Lambda[[r]] = Vartheta[[r]]*Lambda[[r]]
      }
   }
   return(list(Lambda=Lambda, Varphi=Varphi, Vartheta=Vartheta, W=W))
}
