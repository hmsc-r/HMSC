#' @title plotGamma
#'
#' @description Plots heatmaps of parameter estimates or posterior support values of trait effects
#' on species' environmental responses, i.e. how environmental responses in \code{Beta} responds to
#' covariates in \code{X}
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param post posterior summary of Gamma parameters obtained from \code{\link{getPostEstimate}}
#' @param param controls which parameter is plotted, current options include "Mean" for posterior mean
#' estimate, "Support" for the level of statistical support measured by posterior probability for a
#' positive or negative response, and "Sign" to indicate whether the response is positive,
#' negative, or neither of these given the chosen \code{supportLevel}
#' @param trOrder controls the ordering of traits, current options are "Original", and "Vector".
#' If trOrder = "Vector", an ordering vector must be provided (see trVector)
#' @param trVector controls the ordering of traits if trOrder = "Vector". If a subset of traits
#' are listed, only those will be plotted
#' @param covOrder controls the ordering of covariates, current options are "Original" and
#' "Vector". If covOrder = "Vector", an ordering vector must be provided (see covVector)
#' @param covVector controls the ordering of covariates if covOrder = "Vector". If a subset of
#' covariates are listed, only those will be plotted
#' @param trNamesNumbers logical of length 2, where first entry controls whether trait names
#' are added to axes, and second entry controls whether traits numbers are added
#' @param covNamesNumbers logical of length 2, where first entry controls whether covariate names
#' are added to axes, and second entry controls whether covariate numbers are added
#' @param supportLevel controls threshold posterior support for plotting
#' @param cex controls character expansion (font size). Three values, controlling covariate names,
#' trait names, and color legend axis labels
#' @param colors controls the colors of the heatmap, default value \code{colorRampPalette(c("blue","white","red"))}
#' @param colorLevels number of color levels used in the heatmap
#' @param mar plotting margins
#' @param smallplot passed to \code{\link{image.plot}}
#' @param bigplot passed to \code{\link{image.plot}}
#' @param newplot set to  false if the plot will be part of multi-panel plot
#'
#' @examples
#' # Plot posterior support values of trait effects on environmental responses
#' gammaPost=getPostEstimate(TD$m, "Gamma")
#' plotGamma(TD$m, post=gammaPost, param="Support")
#'
#' # Plot parameter estimates of trait effects on environmental responses
#' gammaPost=getPostEstimate(TD$m, "Gamma")
#' plotGamma(TD$m, post=gammaPost, param="Mean")
#'
#' @importFrom graphics par plot.new axis text
#' @importFrom grDevices colorRampPalette
#' @importFrom fields image.plot
#'
#' @export

plotGamma=function(hM, post, param = "Support", trOrder="Original",
  trVector= NULL, covOrder="Original", covVector=NULL, trNamesNumbers=c(TRUE,TRUE),
  covNamesNumbers=c(TRUE,TRUE), supportLevel=.9, cex=c(.8,.8,.8),
  colors=colorRampPalette(c("blue","white","red")), colorLevels = NULL,
  mar=c(6,9,2,0),
  smallplot=NULL, bigplot=NULL, newplot=TRUE){

   if(is.null(colorLevels)){
      if(param=="Sign"){
         colorLevels=3} else {
            colorLevels=200
         }
   }
  switch(class(hM$X)[1L],
         matrix = {
           ncolsX = ncol(hM$X)
         },
         list = {
           ncolsX = ncol(hM$X[[1]])
         }
  )

   covNames = character(hM$nc)
   for (i in 1:hM$nc) {
      sep = ""
      if (covNamesNumbers[1]) {
         covNames[i] = paste(covNames[i], hM$covNames[i], sep = sep)
         sep = " "
      }
      if (covNamesNumbers[2]) {
         covNames[i] = paste(covNames[i], sprintf("(C%d)",
            i), sep = sep)
      }
   }
   trNames = character(hM$nt)
   for (i in 1:hM$nt) {
      sep = ""
      if (trNamesNumbers[1]) {
         trNames[i] = paste(trNames[i], hM$trNames[i], sep = sep)
         sep = " "
      }
      if (trNamesNumbers[2]) {
         trNames[i] = paste(trNames[i], sprintf("(T%d)", i),
            sep = sep)
      }
   }


   if(covOrder=="Vector"){covorder=covVector}
   if(covOrder=="Original"){covorder=1:ncolsX}

   if(trOrder=="Vector"){trorder=trVector}
   if(trOrder=="Original"){trorder=1:ncol(hM$Tr)}


   mgamma=post$mean
   gammaP=post$support

   if(param=="Sign"){
      toPlot = sign(mgamma)
      toPlot = toPlot * ((gammaP>supportLevel) + (gammaP<(1-supportLevel))>0)
      gammaMat = matrix(toPlot, nrow=ncolsX, ncol=ncol(hM$Tr))
   }
   if(param=="Mean"){
      toPlot = mgamma
      toPlot = toPlot * ((gammaP>supportLevel) + (gammaP<(1-supportLevel))>0)
      gammaMat = matrix(toPlot, nrow=ncolsX, ncol=ncol(hM$Tr))
   }
   else{
      if(param=="Support"){
         toPlot = 2*gammaP-1
         toPlot = toPlot * ((gammaP>supportLevel) + (gammaP<(1-supportLevel))>0)
         gammaMat = matrix(toPlot, nrow=ncolsX, ncol=ncol(hM$Tr))
      }}

   rownames(gammaMat) = covNames
   colnames(gammaMat) = trNames
   X = gammaMat[covorder,trorder]

   old.par = par(no.readonly = TRUE)
   colors = colors(colorLevels)

   START=0
   END=.65
   ADJy=1/(ncol(X)*2)
   ADJx=1/(nrow(X)*4)

   if(newplot){
      par(fig = c(0,1,0,1),  mar = mar)
   } else {
      par(old.par, mar=mar)
   }

   plot.new()
    axis(1, at = seq(START+ADJx, END-ADJx,
                    by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)),
         labels = FALSE)
   axis(2, at = seq(ADJy, 1-ADJy, length.out=ncol(X)), labels = FALSE)

   text(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[2],
      labels = covNames[covorder], xpd = TRUE)
   text(y = seq(ADJy, 1-ADJy, length.out=ncol(X)), par("usr")[3] - 0.05, srt = 0, adj = 1,cex = cex[1],
      labels = trNames[trorder], xpd = TRUE)

   if(all(is.na(X)) || sum(abs(X))==0){
      warning("Nothing to plot at this level of posterior support")
      zlim = c(-1,1)
   } else{
      zlim = c(-max(abs(range(X))),max(abs(range(X))))
   }


   image.plot(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)),
      y = seq(ADJy, 1-ADJy, length.out = ncol(X)),
      z = X, add = TRUE, nlevel = colorLevels,
      legend.width = 2, legend.mar = NULL,
      legend.cex = cex,
      axis.args=if(param=="Sign")
      {list(labels=c("+","0","-"),at=c(1,0,-1),cex.axis=cex[3],mgp=c(3,2,0),hadj=1)
      } else {
         list(cex.axis=cex[3],mgp=c(3,2,0),hadj=1)
      },
      graphics.reset = TRUE, horizontal = FALSE, bigplot = bigplot, smallplot = smallplot,
      legend.only = FALSE, col = colors,
      lab.breaks = NULL, zlim = zlim)

   if(newplot){
      par(old.par)
   }
}

