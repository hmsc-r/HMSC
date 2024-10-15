#' @title importPosteriorFromHPC
#'
#' @description Imports the posterior calculated with \code{Hmsc-HPC} add-on to the corresponding Hmsc object
#'
#' @param m the \code{Hmsc} model object used for creating export to Hmsc-HPC fitting
#' @param postList list of MCMC chains, containing the recorded samples from the MCMC sampling
#' @param nSamples the number of recorded MCMC samples per chain, as was used for Hmsc-HPC fitting
#' @param thin the MCMC thinning value, as was used for Hmsc-HPC fitting
#' @param transient the MCMC transient value, as was used for Hmsc-HPC fitting
#' @param alignPost bool flag whether to attempt to align the latent factors and loadings from different chains
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
