#' @title plotBeta
#'
#' @description Plots heatmaps of parameter estimates or posterior support values of species' environmental responses, i.e. how species in \code{Y} responds to covariates in \code{X}
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param post posterior summary of Beta parameters obtained from \code{\link{getPostEstimate}}
#' @param param controls which parameter is plotted, current options include "Mean" for posterior mean estimate,
#' "Support" for the level of statistical support measured by posterior probability for a positive or negative
#' response, and "Sign" to indicate whether the response is positive, negative, or neither of these given the
#' chosen \code{supportLevel}
#' @param plotTree logical. Whether species' environmental responses is to be mapped onto the phylogeny used in
#' model fitting
#' @param SpeciesOrder controls the ordering of species, current options are "Original", "Tree", and "Vector".
#' If SpeciesOrder = "Vector", an ordering vector must be provided (see SpVector). If plotTree = TRUE, SpeciesOrder
#' is ignored
#' @param SpVector controls the ordering of species if SpeciesOrder = "Vector". If a subset of species are listed,
#' only those will be plotted. For alphabetic ordering, try
#' \code{match(1:hM$ns, as.numeric(as.factor(colnames(hM$Y))))}
#' @param covOrder controls the ordering of covariates, current options are "Original" and "Vector".
#' If covOrder = "Vector", an ordering vector must be provided (see covVector)
#' @param covVector controls the ordering of covariates if covOrder = "Vector". If a subset of covariates are listed,
#' only those will be plotted
#' @param spNamesNumbers logical of length 2, where first entry controls whether species names are added to axes,
#' and second entry controls whether species numbers are added
#' @param covNamesNumbers logical of length 2, where first entry controls whether covariate names are added to axes,
#' and second entry controls whether covariate numbers are added
#' @param supportLevel controls threshold posterior support for plotting
#' @param split if plotTree = TRUE, controls the division of the plotting window between the tree and the heatmap.
#' @param cex controls character expansion (font size). Three values, controlling covariate names, species names,
#' and color legend axis labels
#' @param colors controls the colors of the heatmap, default value \code{colorRampPalette(c("blue","white","red"))}
#' @param colorLevels number of color levels used in the heatmap
#' @param mar plotting margins
#' @param marTree plotting margins for phylogenetic tree
#' @param mgp can be used to set the location of the scale bar
#' @param main main title for the plot.
#' @param smallplot passed to \code{\link{image.plot}}
#' @param bigplot passed to \code{\link{image.plot}}
#' @param newplot set to  false if the plot will be part of multi-panel plot initialized with par(mfrow)
#'
#'
#'
#' @examples
#' # Plot posterior support values of species' environmental responses
#' betaPost=getPostEstimate(TD$m, "Beta")
#' plotBeta(TD$m, post=betaPost, param="Support")
#'
#' # Plot parameter estimates of species' environmental responses together with the phylogenetic tree
#' betaPost=getPostEstimate(TD$m, "Beta")
#' plotBeta(TD$m, post=betaPost, param="Mean", plotTree=TRUE)
#'
#' @importFrom graphics par plot plot.new axis text title
#' @importFrom grDevices colorRampPalette
#' @importFrom ape keep.tip read.tree write.tree
#' @importFrom fields image.plot
#' @export


plotBeta = function(hM, post, param = "Support", plotTree = FALSE,
                    SpeciesOrder = "Original", SpVector = NULL,
                    covOrder="Original",
                    covVector=NULL, spNamesNumbers = c(TRUE, TRUE),
                    covNamesNumbers = c(TRUE, TRUE),
                    supportLevel = 0.9, split = 0.3, cex = c(0.7,0.7,0.8),
                    colors = colorRampPalette(c("blue","white","red")), colorLevels = NULL,
                    mar=NULL, marTree=c(6,0,2,0),mgp=c(3,2,0),
                    main = NULL,
                    smallplot=NULL, bigplot=NULL,newplot=TRUE)
{
    ## Check that text arguments are acceptable, and expand to full
    ## names if abbreviated.
    param <- match.arg(param, c("Mean", "Support", "Sign"))
    SpeciesOrder <- match.arg(SpeciesOrder, c("Original", "Vector", "Tree"))
    covOrder <- match.arg(covOrder, c("Original", "Vector"))

   if(is.null(colorLevels)){
      if(param=="Sign"){
         colorLevels=3} else {
            colorLevels=200
         }
   }
   if(is.null(mar)){
      if(sum(spNamesNumbers)>0){
         mar=c(6,6,2,0)
      } else {
         mar=c(6,0,2,0)
      }
   }

   if(plotTree || SpeciesOrder == "Tree"){
      untangle<-function(tree){
         if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
         obj<-attributes(tree)
         tree <- ape::read.tree(text=ape::write.tree(tree))
         ii <- !names(obj) %in% names(attributes(tree))
         attributes(tree)<-c(attributes(tree),obj[ii])
         tree
       }
      tree = keep.tip(hM$phyloTree,hM$spNames)
      tree = untangle(tree)
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

    if(plotTree){
        order=tree$tip.label
    } else {
        order <-
            switch(SpeciesOrder,
                   "Vector" = SpVector,
                   "Original" = rev(1:ncol(hM$Y)),
                   "Tree" = match(tree$tip.label,colnames(hM$Y)))
    }
   if(covOrder=="Vector"){covorder=covVector}
   else if(covOrder=="Original"){covorder=1:hM$nc}

   mbeta=post$mean
   betaP=post$support

   if(param=="Sign"){
      toPlot = sign(mbeta)
      toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
      betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
   }
   else if(param=="Mean"){
      toPlot = mbeta
      toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
      betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
   }
   else{ # param == "Support"
       toPlot = 2*betaP-1
       toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
       betaMat = matrix(toPlot, nrow=hM$nc, ncol=ncol(hM$Y))
   }

   rownames(betaMat) = covNames
   if(plotTree){colnames(betaMat) = gsub(" ", "_", hM$spNames)}
   if(!plotTree){colnames(betaMat) = spNames}

   X = t(betaMat[covorder,order])

   old.par = par(no.readonly = TRUE)
   colors = colors(colorLevels)

   #With tree
   if(plotTree){
      if(newplot){
         par(fig = c(0,split[1],0,1), mar=marTree)
      } else {
         plot.new()
         cpar = par(no.readonly = TRUE)
         cfig = cpar$fig
         nfig = cfig
         nfig[2] = cfig[1]+split[1]*(cfig[2]-cfig[1])
         par(fig=nfig, mar=marTree, new=TRUE)
      }

      if(sum(!spNamesNumbers)==2){plot(tree, show.tip.label=FALSE)}
      else{tree$tip.label[match(gsub(" ", "_", hM$spNames),tree$tip.label)]=spNames
      plot(tree, show.tip.label=TRUE, adj=1, align.tip.label=TRUE, cex=cex[2])}

      if(newplot){
          par(fig = c(split[1],1,0,1),  mar=mar, new=TRUE)
      } else {
         nfig = cfig
         nfig[1] = cfig[1]+split[1]*(cfig[2]-cfig[1])
         par(fig=nfig, new=TRUE)
         par(old.par, mar=mar, new=TRUE)
      }

      START = .05
      END = .7
      ADJy=0
      ADJx=1/(ncol(X)*4)
      plot.new()
      axis(1,seq((START+ADJx), (END-ADJx),
                 by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), labels = FALSE)
      text(x=seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[1],
           labels = covNames[covorder], xpd = TRUE)
      print(par(no.readonly = TRUE)$fig)
   }

   #No tree
   if(!plotTree){
      START=0
      END=.65
      ADJy=1/(nrow(X)*4)
      ADJx=1/(ncol(X)*4)

      plot.new()
      if(newplot){
         par(fig = c(0,1,0,1),  mar = mar)
      } else {
         cpar = par(no.readonly = TRUE)
         cfig = cpar$fig
         par(fig=cfig, mar=mar, new=TRUE)
      }
      axis(1,at = seq((START+ADJx), (END-ADJx),
                      by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)),
           labels = FALSE)
      text(x=seq((START+ADJx), (END-ADJx), by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)), par("usr")[3] - 0.05, srt = 90, adj = 1,cex=cex[1],
           labels = covNames[covorder], xpd = TRUE)
      names=gsub("_"," ",spNames[order])
      text(y = seq(ADJy, 1-ADJy, length.out=nrow(X)),par("usr")[3] - 0.05, srt = 0, adj = 1,cex=cex[2],
           labels = as.expression(lapply(names, function(names) bquote(italic(.(names))))), xpd = TRUE)
   }

   #Plot
   if(all(is.na(X)) || sum(abs(X))==0){
      warning("Nothing to plot at this level of posterior support")
      zlim = c(-1,1)
   } else{
      zlim = c(-max(abs(range(X))),max(abs(range(X))))
   }

   image.plot(x = seq(START+ADJx, END-ADJx, by = ((END-ADJx) - (START+ADJx))/(ncol(X) - 1)),
              y = seq(ADJy, 1-ADJy, length.out=nrow(X)),
              z = t(X), add = TRUE, nlevel = colorLevels, box=TRUE,
              legend.width = 2, legend.mar = NULL,
              legend.cex=cex,
              axis.args=if(param=="Sign")
              {list(labels=c("+","0","-"),at=c(1,0,-1),cex.axis=cex[3],mgp=mgp,hadj=1)
              } else {
                 list(cex.axis=cex[3],mgp=mgp,hadj=1)
              },
              graphics.reset = TRUE, horizontal = FALSE, bigplot = bigplot, smallplot = smallplot,
              legend.only = FALSE, col = colors,zlim = zlim)
   if (!is.null(main))
      title(main = main)

   if(newplot){
      par(old.par)
   } else {
      par(mfrow=cpar$mfrow)
      cmfg=cpar$mfg
      if(cmfg[2]<cmfg[4]){
         cmfg[2] = cmfg[2]+1
      } else {
         if(cmfg[1]<cmfg[3]){
            cmfg[2]=1
            cmfg[1] = cmfg[1]+1
         }
         else{
            cmfg[1]=1
            cmfg[2]=1
         }
      }
      par(mfg=cmfg,mar=cpar$mar)
   }
}

