#' @title coralCombine
#'
#' @description Extracts CORAL-like summary from backbone fitted \code{Hmsc}-class object and combines with CORAL models
#'
#' @param m fitted \code{Hmsc}-class object
#' @param muList.coral matrix of CORAL posterior means
#' @param VList.coral matrix of CORAL posterior flattened variances
#'
#' @return
#' list with combined means and covariance matrices
#'
#'
#' @export

coralCombine = function(m, muList.coral, VList.coral){
   tmp = getPostEstimate(hM=m, parName="Beta")
   mu.backbone = t(tmp$mean)
   colnames(mu.backbone) = m$covNames
   post = poolMcmcChains(m$postList)
   postN = length(post)
   EBeta2.backbone = array(0, dim=c(m$ns,m$nc,m$nc))
   for(i in 1:postN){
      BT = t(post[[i]]$Beta)
      BTA1 = array(BT[,rep(1:m$nc,each=m$nc)], c(m$ns,m$nc,m$nc))
      BTA2 = array(BT[,rep(1:m$nc,m$nc)], c(m$ns,m$nc,m$nc))
      EBeta2.backbone = EBeta2.backbone + BTA1*BTA2
   }
   EBeta2.backbone = EBeta2.backbone / postN
   MA1 = array(mu.backbone[,rep(1:m$nc,each=m$nc)], c(m$ns,m$nc,m$nc))
   MA2 = array(mu.backbone[,rep(1:m$nc,m$nc)], c(m$ns,m$nc,m$nc))
   V.backbone = EBeta2.backbone - MA1*MA2
   mu.all = Reduce(rbind, c(list(mu.backbone), muList.coral))
   V.all = Reduce(rbind, c(list(matrix(V.backbone, m$ns, m$nc^2)), VList.coral))
   return(list(mu=mu.all, V=V.all))
}
