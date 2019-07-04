#' @title biPlot
#'
#' @description Constructs a ordination biplot based on the fitted model
#' @param etaPost Posterior distribution of site loadings Eta
#' @param lambdaPost Posterior distribution of species loadings Lambda
#' @param factors Indices of the two factors to be plotted
#' @param colvar The environmental covariate from XData according to which the sites are to be coloured
#'
#' @examples
#' # Construct an ordination biplot using the first 2 latent factors from a previously fitted HMSC model
#' etaPost = getPostEstimate(TD$m, "Eta")
#' lambdaPost=getPostEstimate(TD$m, "Lambda")
#' biPlot(TD$m, etaPost = etaPost, lambdaPost=lambdaPost, factors=c(1,2))
#'
#'
#' @importFrom graphics plot points text
#' @importFrom grDevices colorRampPalette
#'
#' @export

biPlot=function(hM, etaPost, lambdaPost, factors=c(1,2), colVar=NULL, spNames=hM$spNames, ...){
   if(!is.null(colVar)){
      col = hM$XData[,colVar]
      if (!class(col)=="factor"){
         cols=colorRampPalette(c("blue","white","red"))(hM$np)
         plotorder=order(hM$XData[,which(colnames(hM$XData)==colVar)])
      } else
      {
         cols = col
         plotorder = 1:hM$ny
      }
   }
   else{
      cols="grey"
      plotorder=1:nrow(hM$XData)
   }
   scale1 = abs(c(min(etaPost$mean[,factors[1]]),max(etaPost$mean[,factors[1]])))
   scale2 = abs(c(min(etaPost$mean[,factors[2]]),max(etaPost$mean[,factors[2]])))
   scale1 = min(scale1/abs(c(min(lambdaPost$mean[factors[1],]),max(lambdaPost$mean[factors[1],]))))
   scale2 = min(scale2/abs(c(min(lambdaPost$mean[factors[2],]),max(lambdaPost$mean[factors[2],]))))
   plot(etaPost$mean[,factors[1]][plotorder], etaPost$mean[,factors[2]][plotorder],pch=16, col=cols,
        xlab=paste("Latent variable", factors[1]), ylab=paste("Latent variable", factors[2]), ...)
   points(scale1*lambdaPost$mean[factors[1],], scale2*lambdaPost$mean[factors[2],],pch=17, cex=1)
   text(scale1*lambdaPost$mean[factors[1],], scale2*lambdaPost$mean[factors[2],], spNames, pos=1, cex=1)
}
