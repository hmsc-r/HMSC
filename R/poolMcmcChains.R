#' @title poolMcmcChains
#'
#' @description Combines a list of single or several MCMC chains into a single chain
#'
#' @param postList list of posterior chains
#' @param chainIndex index of chains to be included
#' @param start index of first MCMC sample included
#' @param thin thinning between included MCMC samples
#'
#'
#' @return a list with combined MCMC samples
#'
#' @export

poolMcmcChains = function(postList, chainIndex=1:length(postList), start=1, thin=1){
   post = list()
   for(i in chainIndex){
      chain = postList[[i]]
      ind = seq(start, length(chain), thin)
      post = c(post, chain[ind])
   }
   return(post)
}
