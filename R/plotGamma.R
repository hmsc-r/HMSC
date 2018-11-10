#' @title plotBeta
#'
#' @description Plots betas
#' @param post
#' @param param
#' @param trOrder
#' @param trVector
#' @param covOrder
#' @param covVector
#' @param trNamesNumbers
#' @param covNamesNumbers
#' @param supportLevel
#' @param cex
#'
#' @examples
#'

plotGamma=function(post, param = "Gamma", trOrder="Original", 
  trVector= NULL, covOrder="Original",covVector=NULL, trNamesNumbers=c(T,T),
  covNamesNumbers=c(T,T),supportLevel=.9,cex=c(.8,.8,.8)){
  
   m = self

   covNames = character(m$nc)
   for (i in 1:m$nc) {
      sep = ""
      if (covNamesNumbers[1]) {
         covNames[i] = paste(covNames[i], m$covNames[i], sep = sep)
         sep = " "
      }
      if (covNamesNumbers[2]) {
         covNames[i] = paste(covNames[i], sprintf("(C%d)",
            i), sep = sep)
      }
   }
   trNames = character(m$nt)
   for (i in 1:m$nt) {
      sep = ""
      if (trNamesNumbers[1]) {
         trNames[i] = paste(trNames[i], m$trNames[i], sep = sep)
         sep = " "
      }
      if (trNamesNumbers[2]) {
         trNames[i] = paste(trNames[i], sprintf("(T%d)", i),
            sep = sep)
      }
   }


   if(covOrder=="Vector"){covorder=covVector}
   if(covOrder=="Original"){covorder=1:ncol(m$X)}

   if(trOrder=="Vector"){trorder=trVector}
   if(trOrder=="Original"){trorder=1:ncol(m$Tr)}


   mgamma=post$mean
   gammaP=post$support

   if(param=="Gamma"){
      toPlot = mgamma
      toPlot = toPlot * ((gammaP>supportLevel) + (gammaP<(1-supportLevel))>0)
      gammaMat = matrix(toPlot, nrow=ncol(m$X), ncol=ncol(m$Tr))
   }
   else{
      if(param=="Support"){
         toPlot = 2*gammaP-1
         toPlot = toPlot * ((gammaP>supportLevel) + (gammaP<(1-supportLevel))>0)
         gammaMat = matrix(toPlot, nrow=ncol(m$X), ncol=ncol(m$Tr))
      }}

   rownames(gammaMat) = covNames
   colnames(gammaMat) = trNames
   X = gammaMat[covorder,trorder]

   old.par = par(no.readonly = TRUE)
   colors = colorRampPalette(c("blue","white","red"))(200)

   START=0
   END=.65
   ADJy=1/(ncol(X)*2)
   ADJx=1/(nrow(X)*4)

   par(fig = c(0,1,0,1),  mar = c(6,9,2,0))
   plot.new()
   axis(1,at = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)), labels = F)
   axis(2,at = seq(ADJy, 1-ADJy, length.out=ncol(X)), labels = F)

   text(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[2],
      labels = covNames[covorder], xpd = TRUE)
   text(y = seq(ADJy, 1-ADJy, length.out=ncol(X)), par("usr")[3] - 0.05, srt = 0, adj = 1,cex = cex[1],
      labels = trNames[trorder], xpd = TRUE)

   if(all(is.na(X)) || sum(X)==0){
      warning("Nothing to plot at this level of posterior support")
      zlim = c(-1,1)
   } else{
      zlim = c(-max(abs(range(X))),max(abs(range(X))))
   }


   image.plot(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(nrow(X) - 1)),
      y = seq(ADJy, 1-ADJy, length.out = ncol(X)),
      z = X, add = TRUE, nlevel = 200,
      legend.width = 2, legend.mar = NULL,
      legend.cex = cex, axis.args = list(cex.axis = cex[3], mgp = c(3,2,0), hadj = 1),
      graphics.reset = TRUE, horizontal = FALSE, bigplot = NULL,
      smallplot = NULL, legend.only = FALSE, col = colors,
      lab.breaks = NULL, zlim = zlim)

   par(old.par)
}

Hmsc$set("public", "plotGamma", plotGamma, overwrite=TRUE)

