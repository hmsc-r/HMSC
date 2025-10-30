#' @title importPosteriorFromHPC
#'
#' @description Integrates Hmsc model with a list of chains containing posterior MCMC samples generated with Hmsc-HPC
#'
#' @param m a \code{Hmsc} model object for which the posterior will be imported
#' @param postList list of MCMC chains, where each chain is a list of posterior samples.
#' @param nSamples number of posterior MCMC samples per chain used in Hmsc-HPC
#' @param thin thinning interval used in Hmsc-HPC
#' @param transient number of transient MCMC steps used in Hmsc-HPC
#' @param adaptNf number of MCMC steps for adjustung number of factors used in Hmsc-HPC
#' @param verbose level of verbose to fake
#' @param alignPost boolean flag indicating whether the posterior of each chains should be aligned
#'
#' @return a fitted Hmsc model object
#'
#' @details Naturally, the model object \code{m} must be compatible with the posterior samples imported from Hmsc-HPC.
#' If \code{m} already contains posterior samples, these will be overwritten.
#'
#' Parameters \code{nSamples}, \code{thin}, \code{transient}, \code{adaptNf} and especially \code{verbose} are required
#' to guarantee that the resulted Hmsc object works interchangeably with post-processing functions in the Hmsc package
#' as if the sampling was done in the \code{Hmsc} R package. These values shall mirror the corresponding values
#' that would be used in the \code{sampleMcmc(...)} call.
#'
#' Parameter \code{alignPost} has an identical effect to its counterpart from \code{sampleMcmc(...)}.
#'
#' @examples
#' if(FALSE){
#' # If all chains were run with a single Hmsc-HPC call
#' post_file_path = ... # path to the file where the posterior samples are stored
#' importFromHPC = from_json(readRDS(file = post_file_path)[[1]])
#' postList = importFromHPC[1:nChains]
#' cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))
#' fm = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
#'
#' # If Hmsc-HPC calls were configured to run one chain per call
#' post_file_path = ... # directory with Hmsc-HPC export files stored as post_chainXX_file.rds
#' chainList = vector("list", nChains)
#' for(cInd in 1:nChains){
#'    chain_file_path = file.path(post_file_path, sprintf("post_chain%.2d_file.rds", cInd-1))
#'    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
#' }
#' fm = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
#' }
#'
#' @export


importPosteriorFromHPC = function(m, postList, nSamples, thin, transient, adaptNf=rep(transient,m$nr), verbose=0, alignPost=TRUE){
   m$samples = nSamples
   m$thin = thin
   m$transient = transient
   m$postList = vector("list", length(postList))
   m$adaptNf = adaptNf
   m$verbose = verbose
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
