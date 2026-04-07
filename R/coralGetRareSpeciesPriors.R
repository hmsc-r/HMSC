#' @title coralGetRareSpeciesPriors
#'
#' @description Constructs CORAL priors
#'
#' @param m fitted \code{Hmsc}-class object
#' @param spNames.rare vector of species names that are considered rare
#' @param TrData.rare dataframe of traits for rare species, if non-trivial traits were used in the backbone Hmsc model
#' @param phyloTree phylogeny tree covering both common and rare species
#' @param C.common.rare part of phylogeny similarity matrix between common (columns) and rare (rows) species
#' @param interceptVarInflation multiplier for conditional variance prior of intercept term in the CORAL prior,
#' compared to phylogeny-induced conditional prior
#' @param phyloFast logical, if TRUE use fast tree traversal for phylogeny-induced priors
#'
#' @return
#' A list containing CORAL prior as matrix of means and matrix of flattened covariances
#'
#' @details
#' The returned CORAL prior components are matrices with rows corresponding to the rare species from the \code{spNames.rare} argument.
#'
#' @importFrom ape vcv
#'
#' @export

coralGetRareSpeciesPriors = function(m, spNames.rare, TrData.rare=NULL, phyloTree=NULL, C.common.rare=NULL,
                                     interceptVarInflation=5, phyloFast=FALSE, showProgress=FALSE){
   ns.rare = length(spNames.rare)
   if(is.null(TrData.rare)){
      Tr.rare = matrix(1, ns.rare, 1)
   } else{
      Tr.rare = model.matrix(m$TrFormula, TrData.rare)
   }
   if(phyloFast){
      if(is.null(phyloTree)) phyloTree = m$phyloTree
      spNames.all = c(m$spNames, spNames.rare)
      ns.all = length(spNames.all)
      # phyloTree = keepTipRoot(phyloTree, spNames.all)

      treeList = vector("list", ns.all + phyloTree$Nnode)
      for(i in 1:length(treeList)){
         treeList[[i]] = list(n=0, child=NULL, edgeLen=NULL, parent=0, parentEdgeLen=0)
      }
      for(i in 1:nrow(phyloTree$edge)){
         parentNode = phyloTree$edge[i,1]
         childNodeOrig = phyloTree$edge[i,2]
         if(childNodeOrig <= ns.all){
            childNode = match(phyloTree$tip.label[childNodeOrig], spNames.all)
         } else childNode = childNodeOrig
         nodeObj = treeList[[parentNode]]
         nodeObj$n = nodeObj$n + 1
         nodeObj$child = c(nodeObj$child, childNode)
         nodeObj$edgeLen = c(nodeObj$edgeLen, phyloTree$edge.length[i])
         treeList[[parentNode]] = nodeObj
         nodeObj = treeList[[childNode]]
         nodeObj$parent = parentNode
         nodeObj$parentEdgeLen = phyloTree$edge.length[i]
         treeList[[childNode]] = nodeObj
      }
      phyloTreeRoot = setdiff(1:length(treeList), phyloTree$edge[,2])
   } else {
      if(is.null(C.common.rare)){
         if(is.null(phyloTree)) phyloTree = m$phyloTree
         spNames.all = c(m$spNames, spNames.rare)
         # phyloTree = keepTipRoot(phyloTree, spNames.all)
         C = vcv(phyloTree, model="Brownian", corr=TRUE)
         C12 = C[m$spNames, spNames.rare]
      } else{
         C12 = C.common.rare[m$spNames, spNames.rare]
      }
   }

   E = matrix(0, ns.rare, m$nc)
   E2 = matrix(0, ns.rare, m$nc^2)
   post = poolMcmcChains(m$postList)
   postN = length(post)
   if(showProgress == TRUE){
      pb = txtProgressBar(max=postN, style=3)
   }
   for(i in seq_len(postN)){
      p = post[[i]]
      beta = p$Beta
      gamma = p$Gamma
      V = p$V
      rho = p$rho
      rbeta = beta - tcrossprod(gamma, m$Tr)

      if(phyloFast){
         XTiDX = rep(0, ns.all)
         XTiDS = matrix(0, m$nc, ns.all)
         XTiDX[1:m$ns] = Inf
         XTiDS[, 1:m$ns] = rbeta

         # 1. Scalable conditional mean
         if(rho > 0){
            BetaLambdaCond = recFunScalar_opt(treeList, phyloTreeRoot, rho, XTiDX, XTiDS, sdMult=0)
            muc = t(BetaLambdaCond[, (m$ns+1):ns.all]) + tcrossprod(Tr.rare, gamma)
         } else{
            muc = tcrossprod(Tr.rare, gamma)
         }
         # 2. Scalable mult calculation. Leave it as such for now
         mult = 1 #computeScalableMult(treeList, phyloTreeRoot, rho, m$ns, ns.all)
      } else {
         Q11 = rho*m$C + (1-rho)*diag(m$ns)
         iQ11 = chol2inv(chol(Q11))
         C21iQ11 = t(C12) %*% iQ11
         muc = rho*C21iQ11%*%t(rbeta) + tcrossprod(Tr.rare, gamma)
         mult = 1 - rho^2 * rowSums(C21iQ11 * t(C12))
      }
      E = E + muc
      vV = as.vector(V)
      E2 = E2 + muc[,rep(1:m$nc,each=m$nc)]*muc[,rep(1:m$nc,m$nc)] + matrix(mult,ns.rare,m$nc^2)*matrix(vV,ns.rare,m$nc^2,byrow=TRUE)
      if(showProgress == TRUE){
         setTxtProgressBar(pb, i)
      }
   }
   if(showProgress == TRUE){
      close(pb)
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
