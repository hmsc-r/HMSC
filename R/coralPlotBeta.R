#' @title coralPlotBeta
#'
#' @description Plots summary of Beta coefficents fitted with CORAL approach
#'
#' @param m fitted \code{Hmsc}-class object
#' @param muList.coral arg2
#' @param VList.coral arg3
#'
#' @return
#' list with combined means and covariance matrices
#'
#' @importFrom graphics abline
#' @importFrom grDevices adjustcolor
#'
#' @export

coralPlotBeta = function(mu, V, phyloTree, spNames.common, plotColumns=c(1:ncol(mu)), quantile.support=c(0.05,0.95), plotTree=0.3, seed=NULL,
                         col.common=c("blue","red"), alpha.common=0.5, jitter.common=0.45, cex.common=0.5, pch.common=16,
                         col.rare=c("cyan","pink"), alpha.rare=0.5, jitter.rare=0.2, cex.rare=0.2, pch.rare=16,
                         showCovNames=TRUE){
   if(!is.null(seed)) set.seed(seed)
   spNames = rownames(mu)
   if(!setequal(phyloTree$tip.label, spNames)){
      stop("phyloTree shall contain all species names, provided as rownames(mu), and no other names")
   }

   plotSel = function(k, seln, selp, col, alpha, jitter, cex, pch){
      sel = c(seln, selp)
      co = rep(adjustcolor(col, alpha), c(length(seln), length(selp)))
      rr = sample(length(sel))
      sel = sel[rr]
      co = co[rr]
      points(k+runif(length(sel),-jitter,jitter), sel, cex=cex, pch=pch, col=co)
   }

   ns = length(spNames)
   vars = matrix(NA, ns, ncol(mu))
   V = array(V, c(ns,ncol(mu),ncol(mu)))
   for(i in 1:ncol(mu)){
      vars[,i] = V[,i,i]
   }
   supports = pnorm(q=0, mean=-mu, sd=sqrt(vars))
   rownames(supports) = spNames

   ind = match(phyloTree$tip.label, spNames)
   supportsOrdered = supports[ind,]
   is.common = rownames(supportsOrdered) %in% spNames.common
   is.rare = !is.common

   if(plotTree){
      marTree=c(6,0,2,0)
      old.par = par(no.readonly=TRUE)
      par(fig=c(0.01*plotTree,0.99*plotTree,0.5/(ns+1),(ns+0.5)/(ns+1)), mar=marTree, xaxs="i", yaxs="i")
      plot(phyloTree, show.tip.label=FALSE)
      mar=c(6,0,2,0.2)
      par(fig = c(plotTree,1,0,1), mar=mar, new=TRUE)
   }

   plot(NULL, xlim=c(0.5, length(plotColumns)+0.5), ylim=c(1-0.5,ns+0.5),
        yaxt='n', xaxt='n', ylab=NA, xlab=NA, xaxs ="i", yaxs ="i")
   for(k in 1:length(plotColumns)){
      colInd = plotColumns[k]
      seln.rare = which(is.rare & supportsOrdered[,colInd]<=quantile.support[1])
      selp.rare = which(is.rare & supportsOrdered[,colInd]>=quantile.support[2])
      seln.common = which(is.common & supportsOrdered[,colInd]<=quantile.support[1])
      selp.common = which(is.common & supportsOrdered[,colInd]>=quantile.support[2])
      plotSel(k, seln.rare, selp.rare, col.rare, alpha.rare, jitter.rare, cex.rare, pch.rare)
      plotSel(k, seln.common, selp.common, col.common, alpha.common, jitter.common, cex.common, pch.common)
   }
   abline(v=seq_len(length(plotColumns))+0.5)
   if(showCovNames){
      axis(side=1, at=1:length(plotColumns), labels=colnames(mu)[plotColumns], lwd=0)
   }
   # abline(h=1:ns, lty=2, col="gray")
   # plot(phyloTree, show.tip.label=F)

   if(plotTree){
      par(old.par)
   }
}
# coralPlotBeta(coralAll$mu, coralAll$V, phy, spNames.common, jitter.common=0.01, jitter.rare=0.03, cex.common=1, cex.rare=1)
# coralPlotBeta(t(beta), 0.001*coralAll$V, phy, spNames.common, jitter.common=0.01, jitter.rare=0.03, cex.common=1, cex.rare=1)
