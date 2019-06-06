#' @title plotBeta
#'
#' @description Plots heatmaps of parameter estimates or posterior support values of species' environmental responses, i.e. how species in \code{Y} responds to covariates in \code{X}
#' @param post Posterior summary of Beta parameters obtained from \code{\link{getPostEstimate}}
#' @param param Controls which parameter is plotted, current options include "Mean" for posterior mean estimate, "Support" for the level of statistical support measured by posterior probability for a positive or negative response, and "Sign" to indicate whether the response is positive, negative, or neither of these given the chosen \code{supportLevel}
#' @param plotTree Logical. Whether species' environmental responses is to be mapped onto the phylogeny used in model fitting
#' @param SpeciesOrder Controls the ordering of species, current options are "Original", "Tree", and "Vector". If SpeciesOrder = "Vector", an ordering vector must be provided (see SpVector). If plotTree = T, SpeciesOrder is ignored
#' @param SpVector Controls the ordering of species if SpeciesOrder = "Vector". If a subset of species are listed, only those will be plotted. For alphabetic ordering, try \code{match(1:hM$ns, as.numeric(as.factor(colnames(hM$Y))))}
#' @param covOrder Controls the ordering of covariates, current options are "Original" and "Vector". If covOrder = "Vector", an ordering vector must be provided (see covVector)
#' @param covVector Controls the ordering of covariates if covOrder = "Vector". If a subset of covariates are listed, only those will be plotted
#' @param spNamesNumbers Logical of length 2, where first entry controls whether species names are added to axes, and second entry controls whether species numbers are added
#' @param covNamesNumbers Logical of length 2, where first entry controls whether covariate names are added to axes, and second entry controls whether covariate numbers are added
#' @param supportLevel Controls threshold posterior support for plotting
#' @param split If plotTree = T, controls the division of the plotting window between the tree and the heatmap.
#' @param cex Controls character expansion (font size). Three values, controlling covariate names, species names, and color legend axis labels
#' @param colors Controls the colors of the heatmap, default value \code{colorRampPalette(c("blue","white","red"))}
#'
#'
#' @examples
#' \dontrun{
#' betaPost=getPostEstimate(hM, "Beta")
#' plotBeta(hM, post = betaPost, param = "Support)
#' }
#'
#' @importFrom graphics par plot plot.new axis text
#' @importFrom grDevices colorRampPalette
#' @importFrom ape keep.tip
#' @importFrom fields image.plot
#' @importFrom phytools untangle
#' @export


plotBeta = function(hM, post, param = "Support", plotTree = F,
  SpeciesOrder = "Original", SpVector = NULL, covOrder="Original",
  covVector=NULL, spNamesNumbers = c(T,T), covNamesNumbers = c(T,T),
  supportLevel = 0.9, split = 0.3, cex = c(0.7,0.7,0.8),
  colors = colorRampPalette(c("blue","white","red")),colorLevels = 200,
  smallplot=NULL, bigplot=NULL){

   if(plotTree){
      tree = keep.tip(hM$phyloTree,hM$spNames)
      tree = untangle(tree,"read.tree")
   }

   spNames = character(hM$ns)
   for (i in 1:hM$ns) {
      sep = ""
      if (spNamesNumbers[1]) {
         spNames[i] = paste(spNames[i], hM$spNames[i], sep = sep)
         sep = " "
      }
      if (spNamesNumbers[2]) {
         spNames[i] = paste(spNames[i], sprintf("(S%d)", i),
            sep = sep)
      }
   }

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

   if(plotTree){order=tree$tip.label}
   if(!plotTree & SpeciesOrder=="Vector"){order=SpVector}
   if(!plotTree & SpeciesOrder=="Original"){order=rev(1:ncol(hM$Y))}
   if(!plotTree & SpeciesOrder=="Tree"){order=match(tree$tip.label,colnames(hM$Y))}
   if(!plotTree & SpeciesOrder=="Tree"){order=match(tree$tip.label,hM$spNames)}

   if(covOrder=="Vector"){covorder=covVector}
   if(covOrder=="Original"){covorder=1:hM$nc}

   mbeta=post$mean
   betaP=post$support

   if(param=="Sign"){
      toPlot = sign(mbeta)
      toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
      betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
   }
   if(param=="Mean"){
      toPlot = mbeta
      toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
      betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
   }
   else{
      if(param=="Support"){
         toPlot = 2*betaP-1
         toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
         betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
      }}

   rownames(betaMat) = covNames
   if(plotTree){colnames(betaMat) = gsub(" ", "_", hM$spNames)}
   if(!plotTree){colnames(betaMat) = spNames}

   X = t(betaMat[covorder,order])

   old.par = par(no.readonly = TRUE)
   colors =colors(colorLevels)

   #With tree
   if(plotTree){
      par(fig = c(0,split[1],0,1), mar=c(6,0,2,0))
      if(sum(!spNamesNumbers)==2){plot(tree,show.tip.label=F)}
      else{tree$tip.label[match(gsub(" ", "_", hM$spNames),tree$tip.label)]=spNames
      plot(tree, show.tip.label=T,adj=1,align.tip.label=T,cex=cex[2])}

      par(fig = c(split[1],1,0,1),  mar=c(6,0,2,0), new=T)
      START = .05
      END = .7
      ADJy=0
      ADJx=1/(ncol(X)*4)
      plot.new()
      axis(1,seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), labels = F)
      text(x=seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[1],
         labels = covNames[covorder], xpd = TRUE)
   }

   #No tree
   if(!plotTree){
      START=0
      END=.65
      ADJy=1/(nrow(X)*4)
      ADJx=1/(ncol(X)*4)

      par(fig = c(0,1,0,1),  mar=c(6,10,2,0))
      plot.new()
      axis(1,at = seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), labels = F)
      text(x=seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[1],
         labels = covNames[covorder], xpd = TRUE)
      names=gsub("_"," ",spNames[order])
      text(y = seq(ADJy, 1-ADJy, length.out=nrow(X)),par("usr")[3] - 0.05, srt = 0, adj = 1,cex=cex[2],
         labels = as.expression(lapply(names, function(names) bquote(italic(.(names))))), xpd = TRUE)
   }

   #Plot
   if(all(is.na(X)) || sum(X)==0){
      warning("Nothing to plot at this level of posterior support")
      zlim = c(-1,1)
   } else{
      zlim = c(-max(abs(range(X))),max(abs(range(X))))
   }

   image.plot(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)),
      y = seq(ADJy, 1-ADJy, length.out=nrow(X)),
      z = t(X), add = TRUE, nlevel = colorLevels, box=T,
      legend.width = 2, legend.mar = NULL,
      legend.cex=cex,
      axis.args=if(param=="Sign")
         {list(labels=c("+","0","-"),at=c(1,0,-1),cex.axis=cex[3],mgp=c(3,2,0),hadj=1)
         } else {
            list(cex.axis=cex[3],mgp=c(3,2,0),hadj=1)
         },
      graphics.reset = TRUE, horizontal = FALSE, bigplot = bigplot, smallplot = smallplot,
      legend.only = FALSE, col = colors,
       lab.breaks=NULL, zlim = zlim)
   par(old.par)
}
