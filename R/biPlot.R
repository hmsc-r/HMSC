#' @title biPlot
#'
#' @description Constructs an ordination biplot based on the fitted model
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param etaPost posterior distribution of site loadings (Eta)
#' @param lambdaPost posterior distribution of species loadings
#'     (Lambda)
#' @param factors indices of the two factors to be plotted
#' @param colVar the environmental covariate from XData according to
#'     which the sites are to be coloured
#' @param colors controls the colors of the heatmap. For continuous
#'     covariates, colors should be given as a name of palette, with
#'     default value \code{colorRampPalette(c("blue","white","red"))},
#'     or as a vector of colours. For factors, colors should be given
#'     as a vector of colours, e.g. \code{c("blue","red")}.
#' @param spNames a vector of species names to be added to the
#'     ordination diagram
#' @param \dots other parameters passed to the function.
#'
#' @examples
#' # Construct an ordination biplot using two chosen latent factors from a previously fitted HMSC model
#' etaPost = getPostEstimate(TD$m, "Eta")
#' lambdaPost=getPostEstimate(TD$m, "Lambda")
#' biPlot(TD$m, etaPost=etaPost, lambdaPost=lambdaPost, factors=c(1,2))
#'
#'
#' @importFrom graphics plot points text
#' @importFrom grDevices colorRampPalette palette
#'
#' @export

biPlot=function(hM, etaPost, lambdaPost, factors=c(1,2), colVar=NULL, colors=NULL, spNames=hM$spNames, ...){
   if(!is.null(colVar)){
      if (!is.null(hM$XData))
          col = hM$XData[,colVar]
      else
          col = hM$X[, colVar]
      if (is.null(col))
          stop("colVar was not found")
      if (!is.factor(col)){
         if(is.null(colors)){
            colors = colorRampPalette(c("blue","white","red"))
         }
         ## evaluate colors to colour vector if a function
         if (is.function(colors))
             colors <- colors(100)
         ## scale colVar to colors index
         col = round((col - min(col))/diff(range(col)) *
                     (length(colors) - 1) + 1)
         cols = colors[col]
      } else {
         if(is.null(colors)){
            colors = palette()
            cols = colors[col]
         }
         cols = colors[col]
      }
   }
   else{
      cols="grey"
   }
   scale1 = abs(c(min(etaPost$mean[,factors[1]]),max(etaPost$mean[,factors[1]])))
   scale2 = abs(c(min(etaPost$mean[,factors[2]]),max(etaPost$mean[,factors[2]])))
   scale1 = min(scale1/abs(c(min(lambdaPost$mean[factors[1],]),max(lambdaPost$mean[factors[1],]))))
   scale2 = min(scale2/abs(c(min(lambdaPost$mean[factors[2],]),max(lambdaPost$mean[factors[2],]))))
   scale <- min(scale1, scale2)
   plot(etaPost$mean[,factors[1]], etaPost$mean[,factors[2]], pch=16, col=cols,
        xlab=paste("Latent variable", factors[1]), ylab=paste("Latent variable", factors[2]),
        asp = 1, ...)
   points(scale*lambdaPost$mean[factors[1],], scale*lambdaPost$mean[factors[2],],pch=17, cex=1)
   text(scale*lambdaPost$mean[factors[1],], scale*lambdaPost$mean[factors[2],], spNames, pos=1, cex=1)
}
