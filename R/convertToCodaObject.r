#' @title Hmsc$convertToCodaObject
#'
#' @description Converts Hmsc posterior structure into a named list of mcmc.list objects
#' @param start
#'
#' @examples
#'

convertToCodaObject = function(start=1, spNamesNumbers=c(TRUE,TRUE), covNamesNumbers=c(TRUE,TRUE), trNamesNumbers=c(TRUE,TRUE),
   Beta=TRUE, Gamma=TRUE, V=TRUE, Sigma=TRUE, Rho=TRUE, Eta=TRUE, Lambda=TRUE, Alpha=TRUE,
   Omega=TRUE, Psi=TRUE, Delta=TRUE) {
   m = self
   if (is.null(m$C)){
      Rho = FALSE
   }

   nChains = length(m$postList)
   nc = m$nc
   nt = m$nt
   ns = m$ns
   nr = m$nr

   thin = m$thin
   samples = m$samples
   start1 = m$transient+start*thin
   end1 = m$transient + samples*thin


   spNames = character(ns)
   for (i in 1:ns){
      sep = ""
      if (spNamesNumbers[1]){
         spNames[i] = paste(spNames[i],m$spNames[i],sep=sep)
         sep = " "
      }
      if (spNamesNumbers[2]){
         spNames[i] = paste(spNames[i],sprintf("(S%d)",i),sep=sep)
      }
   }

   trNames = character(nt)
   for (i in 1:nt){
      sep = ""
      if (trNamesNumbers[1]){
         trNames[i] = paste(trNames[i],m$trNames[i],sep=sep)
         sep = " "
      }
      if (trNamesNumbers[2]){
         trNames[i] = paste(trNames[i],sprintf("(T%d)",i),sep=sep)
      }
   }

   covNames = character(nc)
   for (i in 1:nc){
      sep = ""
      if (covNamesNumbers[1]){
         covNames[i] = paste(covNames[i],m$covNames[i],sep=sep)
         sep = " "
      }
      if (covNamesNumbers[2]){
         covNames[i] = paste(covNames[i],sprintf("(C%d)",i),sep=sep)
      }
   }

   postListAll = m$postList

   postBeta = list()
   postGamma = list()
   postV = list()
   postSigma = list()
   postRho = list()
   postEta = list()
   postLambda = list()
   postOmega = list()
   postAlpha = list()
   postPsi = list()
   postDelta = list()

   for (chain in 1:nChains){
      postList = postListAll[[chain]][start:length(postListAll[[chain]])]

      if (Beta){
         tmp = do.call(rbind, lapply(postList, function(a) as.vector(a$Beta) ))
         colnames(tmp) = sprintf("B[%s, %s]",covNames[rep(1:nc,ns)],spNames[rep(1:ns,each=nc)])
         postBeta[[chain]] = mcmc(tmp, thin=thin, start=start1, end=end1)
      }

      if (Gamma){
         tmp = do.call(rbind, lapply(postList, function(a) as.vector(a$Gamma) ))
         colnames(tmp) = sprintf("G[%s, %s]",covNames[rep(1:nc,nt)],trNames[rep(1:nt,each=nc)])
         postGamma[[chain]] = mcmc(tmp, thin=thin, start=start1, end=end1)
      }

      if (V){
         tmp = do.call(rbind, lapply(postList, function(a) as.vector(a$V) ))
         colnames(tmp) = sprintf("V[%s, %s]",covNames[rep(1:nc,nc)],covNames[rep(1:nc,each=nc)])
         postV[[chain]] = mcmc(tmp, thin=thin, start=start1, end=end1)
      }

      if (Sigma){
         tmp = do.call(rbind, lapply(postList, function(a) a$sigma ))
         colnames(tmp) = sprintf("Sig[%s]",spNames)
         postSigma[[chain]] = mcmc(tmp, thin=thin, start=start1, end=end1)
      }

      if (Rho){
         postRho[[chain]] = mcmc(unlist(lapply(postList, function(a) a$rho)), thin=thin, start=start1, end=end1)
      }

      getOmega = function(a,r=1)
         return(crossprod(a$Lambda[[r]]))
      getEta = function(a,r=1)
         return(a$Eta[[r]])
      getLambda = function(a,r=1)
         return(a$Lambda[[r]])
      getAlpha = function(a,r=1)
         return(a$Alpha[[r]])
      getPsi = function(a,r=1)
         return(a$Psi[[r]])
      getDelta = function(a,r=1)
         return(a$Delta[[r]])

      postEta1 = vector("list", nr)
      postLambda1 = vector("list", nr)
      postOmega1 = vector("list", nr)
      postAlpha1 = vector("list", nr)
      postPsi1 = vector("list", nr)
      postDelta1 = vector("list", nr)

      for(r in 1:nr){
         postNf = unlist(lapply(lapply(postList, getEta, r=r), ncol))
         if(length(unique(postNf))!=1)
            stop("HMSC: number of latent factors was changing in selected sequence of samples")
         nf = unique(postNf)
         if (Eta)      {
            tmp = do.call(rbind, lapply(lapply(postList, getEta, r=r), as.vector))
            colnames(tmp) = sprintf("Eta%d[%s, factor%s]",r,m$dfPi[(rep(1:m$np[r],nf)),r],as.character(rep(1:nf,each=m$np[r])))
            postEta1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Lambda)      {
            tmp = do.call(rbind, lapply(lapply(postList, getLambda, r=r), as.vector))
            colnames(tmp) = sprintf("Lambda%d[%s, factor%s]",r,spNames[rep(1:ns,nf)],as.character(rep(1:nf,each=ns)) )
            postLambda1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }
         if (Omega)      {
            tmp = do.call(rbind, lapply(lapply(postList, getOmega, r=r), as.vector))
            colnames(tmp) = sprintf("Omega%d[%s, %s]",r,spNames[rep(1:ns,ns)],spNames[(rep(1:ns,each=ns))] )
            postOmega1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }
         if (Alpha)      {
            tmp = do.call(rbind, lapply(postList, getAlpha, r=r))
            colnames(tmp) = sprintf("Alpha%d[factor%s]",r,as.character(1:nf) )
            postAlpha1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Psi)    {
            tmp = do.call(rbind, lapply(lapply(postList, getPsi, r=r), as.vector))
            colnames(tmp) = sprintf("Psi%d[%s, factor%s]",r,spNames[rep(1:ns,nf)],as.character(rep(1:nf,each=ns)) )
            postPsi1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Delta)      {
            tmp = do.call(rbind, lapply(lapply(postList, getDelta, r=r), as.vector))
            colnames(tmp) = sprintf("Delta%d[factor%s]",r,as.character(1:nf) )
            postDelta1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }
      }

      postEta[[chain]] = postEta1
      postLambda[[chain]] = postLambda1
      postOmega[[chain]] = postOmega1
      postAlpha[[chain]] = postAlpha1
      postPsi[[chain]] = postPsi1
      postDelta[[chain]] = postDelta1
   }

   mpost=list()

   if (Beta){
      mpost$Beta = mcmc.list(postBeta)
   }
   if (Gamma){
      mpost$Gamma = mcmc.list(postGamma)
   }
   if (V){
      mpost$V = mcmc.list(postV)
   }
   if (Sigma){
      mpost$Sigma = mcmc.list(postSigma)
   }
   if (Rho){
      mpost$Rho = mcmc.list(postRho)
   }

   for(r in 1:nr){
      postEta2 = vector("list", nChains)
      postLambda2 = vector("list", nChains)
      postOmega2 = vector("list", nChains)
      postAlpha2 = vector("list", nChains)
      postPsi2 = vector("list", nChains)
      postDelta2 = vector("list", nChains)

      for (chain in 1:nChains){
         postEta1=postEta[[chain]]
         postEta2[[chain]] = postEta1[[r]]

         postLambda1=postLambda[[chain]]
         postLambda2[[chain]] = postLambda1[[r]]

         postOmega1=postOmega[[chain]]
         postOmega2[[chain]] = postOmega1[[r]]

         postAlpha1=postAlpha[[chain]]
         postAlpha2[[chain]] = postAlpha1[[r]]

         postPsi1=postPsi[[chain]]
         postPsi2[[chain]] = postPsi1[[r]]

         postDelta1=postDelta[[chain]]
         postDelta2[[chain]] = postDelta1[[r]]
      }

      if (Eta){
         mpost$Eta[[r]] = mcmc.list(postEta2)
      }
      if (Lambda){
         mpost$Lambda[[r]] = mcmc.list(postLambda2)
      }
      if (Omega){
         mpost$Omega[[r]] = mcmc.list(postOmega2)
      }
      if (Alpha){
         mpost$Alpha[[r]] = mcmc.list(postAlpha2)
      }
      if (Psi){
         mpost$Psi[[r]] = mcmc.list(postPsi2)
      }
      if (Delta){
         mpost$Delta[[r]] = mcmc.list(postDelta2)
      }
   }
   return(mpost)
}

Hmsc$set("public", "convertToCodaObject", convertToCodaObject, overwrite=TRUE)


