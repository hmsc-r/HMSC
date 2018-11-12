#' @title poolMcmcChains
#'
#' @description Combines a list MCMC chains into single chain and
#'
#' @param postList list of posterior chains
#' @param chainIndex index of chains to be included
#' @param start index of first MCMC sample included
#' @param thin thinning between included MCMC samples
#'
#'
#' @return
#'
#'
#' @seealso
#'
#' 
#' @examples
#'
#' @export

poolMcmcChains = function(postList, chainIndex=1:length(postList), 
  start=1, thin=1){
   post = list()
   for(i in chainIndex){
      chain = postList[[i]]
      ind = seq(start, length(chain), thin)
      post = c(post, chain[ind])
   }
   return(post)
}
