#' @title biPlot
#'
#' @description Makes a orindation biPlot based on the fitted model
#' @param etaPost Posterior distribution of site loadings Eta
#' @param lambdaPost Posterior distribution of species loadings Lambda
#' @param factors the indices of the two factors to be shown
#' @param colvar the environmental covariate from XData according to which the sites are to be coloured
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' biPlot(m, etaPost = getPostEstimate(m, "Eta"), lambdaPost=getPostEstimate(m, "Lambda"), factors=c(1,2))
#'
#' @export

biPlot=function(hM, etaPost, lambdaPost, factors=c(1,2), colVar=NULL){
   if(!is.null(colVar)){
      cols=colorRampPalette(c("blue","white","red"))(100)
      plotorder=order(hM$XData[,which(colnames(m$XData)==colVar)])
   }
   else{
      cols="grey"
      plotorder=1:nrow(hM$XData)
   }
   plot(etaPost$mean[,factors[1]][plotorder], etaPost$mean[,factors[2]][plotorder],pch=16, col=cols,
        xlab=paste("Latent variable", factors[1]), ylab=paste("Latent variable", factors[2]))
   points(lambdaPost$mean[factors[1],], lambdaPost$mean[factors[2],],pch=17, cex=1)
   text(lambdaPost$mean[factors[1],], lambdaPost$mean[factors[2],], m$spNames, pos=1, cex=1)
}
