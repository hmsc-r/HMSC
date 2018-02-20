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

sampleMcmc = function(samples, thin=1, initPar=NULL, repN=1, saveToDisk=FALSE, verbose=samples*thin/100){
   self$samples = samples
   self$thin = thin
   self$repN = repN
   self$saveToDisk = saveToDisk

   X = self$X
   Tr = self$Tr

   mGamma = self$mGamma
   iUGamma = chol2inv(solve(self$UGamma))
   V0 = self$V0
   f0 = self$f0


   parList = private$computeInitialParameters(initPar)

   Gamma = parList$Gamma
   V = parList$V
   iV = solve(V)
   Beta = parList$Beta
   sigma = parList$sigma
   Z = parList$Z


   repList = vector("list", repN)
   for(repN in 1:repN){
      postList = vector("list", samples)
      for(iter in 1:(samples*thin)){
         Beta = updateBeta(Z,Gamma,V,sigma, X,Tr)
         GammaVList = updateGammaV(Beta,Gamma,iV, Tr, mGamma,iUGamma,V0,f0)
         Gamma = GammaVList$Gamma
         iV = GammaVList$iV
         Z = updateZ(Z, Beta, X)

         if(iter %% thin == 0){
            postList[[iter/thin]] = private$combineParameters(Beta, Gamma, V)
         }
         if(iter %% verbose == 0){
            print(sprintf("Replicate %d, iteration %d of %d", repN, iter, samples*thin))
         }
      }
      repList[[repN]] = postList
      self$postList = postList
   }
   self$repList = repList
}

Hmsc$set("public", "sampleMcmc", sampleMcmc, overwrite=TRUE)

