#' @title sampleMCMC
#'
#' @description Samples the posterior with block-conditional Gibbs MCMC sampler
#'
#' @param samples number of MCMC steps to be recorded
#' @param thin thinning between recorded MCMC samples
#' @param initPar initial parameters value
#' @param replicates number of replicates the samples should be splitted
#' @param saveToDisk whether to save replicates to the disk once they are ready and discard them from RAM memory
#' @param discardDataPar whether to discard the precomputed data-based quantities once the sampling is finished
#'
#' @examples
#'
#' @export

sampleMCMC = function(samples, thin=1, initPar=NULL, repN=1, saveToDisk=FALSE){
   self$samples = samples
   self$thin = thin
   self$repN = repN
   self$saveToDisk = saveToDisk

   parList = self$computeInitialParameters(initPar)

   repList = vector("list", repN)
   for(repN in 1:repN){
      postList = vector("list", samples)
      for(iter in 1:(samples*thin)){
         Beta = updateBeta(Z, Gamma, V, Tr)
         GammaVList = updateGamma(Beta, mGamma, UGamma)
         Gamma = GammaVList$Gamma
         V = GammaVList$V
         Z = updateZ(Z, Beta, X)

         if(iter %% thin == 0){
            parList[[iter/thin]] = self$combinePar(Beta, Gamma, V)
         }
      }
   }
}

Hmsc$set("public", "sample", sampleMCMC, overwrite=TRUE)

