# Developed by Oystein
plotBeta=function(post, param = "Support", plotTree = F, SpeciesOrder = "Original", SpVector = NULL, covOrder="Original",covVector=NULL, spNamesNumbers = c(T,T), covNamesNumbers = c(T,T), supportLevel = .9, split = .3, cex = c(0.7,0.7,0.8)){
   m = self

   if(plotTree){
      tree=m$phyloTree
      tree=untangle(tree,"read.tree")
   }

   spNames = character(m$ns)
   for (i in 1:m$ns) {
      sep = ""
      if (spNamesNumbers[1]) {
         spNames[i] = paste(spNames[i], m$spNames[i], sep = sep)
         sep = " "
      }
      if (spNamesNumbers[2]) {
         spNames[i] = paste(spNames[i], sprintf("(S%d)", i),
            sep = sep)
      }
   }

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

   if(plotTree){order=tree$tip.label}
   if(!plotTree & SpeciesOrder=="Vector"){order=SpVector}
   if(!plotTree & SpeciesOrder=="Original"){order=rev(1:ncol(m$Y))}
   if(!plotTree & SpeciesOrder=="Tree"){order=match(tree$tip.label,colnames(m$Y))}
   if(!plotTree & SpeciesOrder=="Tree"){order=match(tree$tip.label,m$spNames)}

   if(covOrder=="Vector"){covorder=covVector}
   if(covOrder=="Original"){covorder=1:ncol(m$X)}

   mbeta=post$mean
   betaP=post$support

   if(param=="Beta"){
      toPlot = mbeta
      toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
      betaMat = matrix(toPlot, nrow=ncol(m$X), ncol=ncol(m$Y))
   }
   else{
      if(param=="Support"){
         toPlot = 2*betaP-1
         toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
         betaMat = matrix(toPlot, nrow=ncol(m$X), ncol=ncol(m$Y))
      }}

   rownames(betaMat) = covNames
   if(plotTree){colnames(betaMat) = m$spNames}
   if(!plotTree){colnames(betaMat) = spNames}

   X = t(betaMat[covorder,order])

   old.par = par(no.readonly = TRUE)
   colors = colorRampPalette(c("blue","white","red"))(200)

   #With tree
   if(plotTree){
      par(fig = c(0,split[1],0,1), mar=c(6,0,2,0))
      if(sum(!spNamesNumbers)==2){plot(tree,show.tip.label=F)}
      else{tree$tip.label[match(m$spNames,tree$tip.label)]=spNames
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
      z = t(X), add = TRUE, nlevel = 200, box=T,
      legend.width = 2, legend.mar = NULL,
      legend.cex=cex, axis.args=list(cex.axis=cex[3],mgp=c(3,2,0),hadj=1),
      graphics.reset = TRUE, horizontal = FALSE, bigplot = NULL,
      smallplot = NULL, legend.only = FALSE, col = colors,
      lab.breaks=NULL, zlim = zlim)
   par(old.par)
}

Hmsc$set("public", "plotBeta", plotBeta, overwrite=TRUE)
