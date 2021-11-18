#' @title convertToCodaObject
#'
#' @description Converts the Hmsc posterior into a named list of mcmc.list objects
#' @param hM a fitted \code{Hmsc} model object
#' @param start index of first MCMC sample included
#' @param spNamesNumbers logical of length 2, where first entry controls whether species names are printed,
#' and second entry controls whether species numbers are printed
#' @param covNamesNumbers Logical of length 2, where first entry controls whether covariate names are printed,
#' and second entry controls whether covariate numbers are printed
#' @param trNamesNumbers Logical of length 2, where first entry controls whether trait names are printed,
#' and second entry controls whether traits numbers are printed
#' @param Beta logical indicating whether posterior of Beta is included
#' @param Gamma logical indicating whether posterior of Gamma is included
#' @param V logical indicating whether posterior of V is included
#' @param Sigma logical indicating whether posterior of Sigma is included
#' @param Rho logical indicating whether posterior of Rho is included
#' @param Eta logical indicating whether posterior of Eta is included
#' @param Lambda logical indicating whether posterior of Lambda is included
#' @param Alpha logical indicating whether posterior of Alpha is included
#' @param Omega logical indicating whether posterior of Omega is included
#' @param Psi logical indicating whether posterior of Psi is included
#' @param Delta logical indicating whether posterior of Delta is included
#'
#' @return A named list that can be analysed with \CRANpkg{coda} functions.
#'
#' @examples
#' # Convert recorded posterior samples in \code{Hmsc} object to coda object
#' codaObject = convertToCodaObject(TD$m)
#'
#' # Convert recorded posterior samples, starting from sample 100, in m object to coda object
#' codaObject = convertToCodaObject(TD$m, start=100)
#'
#' @importFrom coda mcmc mcmc.list
#' @export

convertToCodaObject = function(hM, start=1, spNamesNumbers=c(TRUE,TRUE),
  covNamesNumbers=c(TRUE,TRUE), trNamesNumbers=c(TRUE,TRUE),
  Beta=TRUE, Gamma=TRUE, V=TRUE, Sigma=TRUE, Rho=TRUE, Eta=TRUE, Lambda=TRUE,
  Alpha=TRUE, Omega=TRUE, Psi=TRUE, Delta=TRUE)
{
   if (is.null(hM$postList))
      stop("Hmsc object ", sQuote(substitute(hM)), " has no posterior samples")
   if (is.null(hM$C)){
      Rho = FALSE
   }

   nChains = length(hM$postList)
   nc = hM$nc
   nt = hM$nt
   ns = hM$ns
   nr = hM$nr

   thin = hM$thin
   samples = hM$samples
   start1 = hM$transient+start*thin
   end1 = hM$transient + samples*thin


   spNames = character(ns)
   for (i in 1:ns){
      sep = ""
      if (spNamesNumbers[1]){
         spNames[i] = paste(spNames[i],hM$spNames[i],sep=sep)
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
         trNames[i] = paste(trNames[i],hM$trNames[i],sep=sep)
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
         covNames[i] = paste(covNames[i],hM$covNames[i],sep=sep)
         sep = " "
      }
      if (covNamesNumbers[2]){
         covNames[i] = paste(covNames[i],sprintf("(C%d)",i),sep=sep)
      }
   }

   postListAll = hM$postList

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

   nfMat = matrix(0,nChains,nr)
   for (chain in 1:nChains){
      postList = postListAll[[chain]][start:length(postListAll[[chain]])]
      for(r in seq_len(nr)){
         nfMat[chain,r] = max(unlist(lapply(lapply(postList, getEta, r=r), ncol)))
      }
   }
   nfMax = apply(nfMat,2,max)

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

      postEta1 = vector("list", nr)
      postLambda1 = vector("list", nr)
      postOmega1 = vector("list", nr)
      postAlpha1 = vector("list", nr)
      postPsi1 = vector("list", nr)
      postDelta1 = vector("list", nr)

      for(r in seq_len(nr)){
         postNf = unlist(lapply(lapply(postList, getEta, r=r), ncol))
         if(length(unique(postNf))!=1)
            stop("HMSC: number of latent factors was changing in selected sequence of samples")
         nf = unique(postNf)
         if (Eta)      {
            tmp1 = lapply(postList, getEta, r=r)
            tmp2 = lapply(tmp1, function(A) cbind(A, matrix(0,nrow(A),nfMax[r]-ncol(A))))
            tmp = do.call(rbind, lapply(tmp2, as.vector))
            colnames(tmp) = sprintf("Eta%d[%s, factor%s]",r,levels(hM$dfPi[,r])[(rep(1:hM$np[r],nfMax[r]))],as.character(rep(1:nfMax[r],each=hM$np[r])))
            postEta1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Lambda)      {
            tmp1 = lapply(postList, getLambda, r=r)
            tmp2 = lapply(tmp1, function(A) rbind(A, matrix(0,nfMax[r]-nrow(A),ncol(A))))
            tmp = do.call(rbind, lapply(tmp2, function(A) as.vector(t(A)) ))
            colnames(tmp) = sprintf("Lambda%d[%s, factor%s]",r,spNames[rep(1:ns,nfMax[r])],as.character(rep(1:nfMax[r],each=ns)) )
            postLambda1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }
         if (Omega)      {
            tmp = do.call(rbind, lapply(lapply(postList, getOmega, r=r), as.vector))
            colnames(tmp) = sprintf("Omega%d[%s, %s]",r,spNames[rep(1:ns,ns)],spNames[(rep(1:ns,each=ns))] )
            postOmega1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }
         if (Alpha)      {
            tmp1 = lapply(lapply(postList, getAlpha, r=r), function(a) hM$rL[[r]]$alphapw[a,1])
            tmp2 = lapply(tmp1, function(a) c(a,rep(0,nfMax[r]-length(a))))
            tmp = do.call(rbind, tmp2)
            colnames(tmp) = sprintf("Alpha%d[factor%s]",r,as.character(1:nfMax[r]) )
            postAlpha1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Psi)    {
            tmp1 = lapply(postList, getPsi, r=r)
            tmp2 = lapply(tmp1, function(A) rbind(A, matrix(0,nfMax[r]-nrow(A),ncol(A))))
            tmp = do.call(rbind, lapply(tmp2, function(A) as.vector(t(A)) ))
            colnames(tmp) = sprintf("Psi%d[%s, factor%s]",r,spNames[rep(1:ns,nfMax[r])],as.character(rep(1:nfMax[r],each=ns)) )
            postPsi1[[r]] = mcmc(tmp, thin=thin, start=start1, end=end1)
         }

         if (Delta)      {
            tmp1 = lapply(lapply(postList, getDelta, r=r), function(a) c(a,rep(0,nfMax[r]-length(a))))
            tmp = do.call(rbind, tmp1)
            colnames(tmp) = sprintf("Delta%d[factor%s]",r,as.character(1:nfMax[r]) )
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

   for(r in seq_len(nr)){
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


