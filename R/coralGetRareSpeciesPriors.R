#' @title coralGetRareSpeciesPriors
#'
#' @description Constructs CORAL priors
#'
#' @param m fitted \code{Hmsc}-class object
#' @param spNames.rare arg2
#' @param TrData.rare arg3
#' @param phyloTree arg4
#' @param C.common.rare arg5
#' @param interceptVarInflation arg6
#'
#' @return
#' A list containing means and covariances of CORAL prior
#'
#' @importFrom ape vcv
#'
#' @export

coralGetRareSpeciesPriors = function(m, spNames.rare, TrData.rare=NULL, phyloTree=NULL, C.common.rare=NULL, interceptVarInflation=5){
   # requested input - fitted Hmsc backbone model and names of rare species
   # additionally if Hmsc model was using traits, a corresponding TrData.rare must be given for rare species
   # finally, if phylogeny given to Hmsc model does not cover the rare species, a joint phylogeny shall be provided

   #TODO We shall avoid construction of C explicitely, as it is HUGE and SLOW once the number of rare is large.
   # Best will be tweak the ape/phytools functionality somehow, but I am unaware of how that can be done.
   # Alternatively, we need to disentangle the phylogeny tree structure manually to get these common-rare correlations.
   # Can be done by searching the tip-closest common ancestor for each common-rare pair,
   # but perhaps some vectorised solution is possible.

   ns.rare = length(spNames.rare)
   if(is.null(TrData.rare)){
      Tr.rare = matrix(1, ns.rare, 1)
   } else{
      Tr.rare = model.matrix(m$TrFormula, TrData.rare)
   }
   if(is.null(C.common.rare)){
      if(is.null(phyloTree)) phyloTree = m$phyloTree
      spNames.all = c(m$spNames, spNames.rare)
      phyloTree = keepTipRoot(phyloTree, spNames.all)
      C = vcv(phyloTree, model="Brownian", corr=TRUE)
      C12 = C[m$spNames, spNames.rare]
   } else{
      C12 = C.common.rare[m$spNames, spNames.rare]
   }
   E = matrix(0, ns.rare, m$nc)
   E2 = matrix(0, ns.rare, m$nc^2)
   post = poolMcmcChains(m$postList)
   postN = length(post)
   for(i in 1:postN){
      p = post[[i]]
      beta = p$Beta
      gamma = p$Gamma
      V = p$V
      rho = p$rho
      Q11 = rho*m$C + (1-rho)*diag(m$ns)
      iQ11 = chol2inv(chol(Q11))
      rbeta = beta - tcrossprod(gamma, m$Tr)
      C21iQ11 = t(C12) %*% iQ11
      muc = rho*C21iQ11%*%t(rbeta) + tcrossprod(Tr.rare, gamma)
      mult = 1 - rho^2 * rowSums(C21iQ11 * t(C12))
      E = E + muc
      vV = as.vector(V)
      E2 = E2 + muc[,rep(1:m$nc,each=m$nc)]*muc[,rep(1:m$nc,m$nc)] + matrix(mult,ns.rare,m$nc^2)*matrix(vV,ns.rare,m$nc^2,byrow=TRUE)
   }
   E = E/postN
   E2 = E2/postN
   mu.coral = E
   V.coral = E2 - E[,rep(1:m$nc,each=m$nc)]*E[,rep(1:m$nc,m$nc)]
   indIntercept = which(m$covNames %in% c("intercept", "Intercept", "(Intercept)"))
   if(length(indIntercept) > 0){
      VArray = array(V.coral, c(ns.rare,m$nc,m$nc))
      VArray[,,indIntercept] = sqrt(interceptVarInflation) * VArray[,,indIntercept]
      VArray[,indIntercept,] = sqrt(interceptVarInflation) * VArray[,indIntercept,]
      V.coral = matrix(VArray, ns.rare, m$nc^2)
   }
   rownames(mu.coral) = rownames(V.coral) = spNames.rare
   return(list(mu=mu.coral, V=V.coral))
}
