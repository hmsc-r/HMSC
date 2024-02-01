#' @title importPosteriorFromHPC
#'
#' @description Computes initial parameter values before the sampling starts
#'
#' @param m a \code{Hmsc} model object
#' @param initPar a list of initial parameter values
#' @export


importPosteriorFromHPC = function(m, postList, nSamples, thin, transient, alignPost=TRUE){
   m$samples = nSamples
   m$thin = thin
   m$transient = transient
   m$postList = vector("list", length(postList))
   for(cInd in 1:length(postList)){
      if(m$samples != length(postList[[cInd]])){
         stop("Hmsc:importPosteriorFromHPC = each chain length miust be equalt to nSamples")
      }
      m$postList[[cInd]] = vector("list", m$samples)
      for(sInd in 1:length(postList[[cInd]])){
         s = postList[[cInd]][[sInd]]
         m$postList[[cInd]][[sInd]] =
            combineParameters(s$Beta, s$BetaSel, s$wRRR, s$Gamma, s$iV, s$rho, s$sigma**-2,
                              s$Eta, s$Lambda, s$Alpha, s$Psi, s$Delta, s$PsiRRR, s$DeltaRRR,
                              m$ncNRRR, m$ncRRR, m$ncsel, m$XSelect, m$XScalePar, m$XInterceptInd,
                              m$XRRRScalePar, m$nt, m$TrScalePar, m$TrInterceptInd, m$rhopw)
      }
   }
   if (alignPost){
      for (i in 1:5){
         m = alignPosterior(m)
      }
   }
   return(m)
}
