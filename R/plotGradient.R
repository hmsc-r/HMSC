#' @title Hmsc$plotGradient
#'
#' @description Plots an environmental gradient over one of the variables included in \code{XData}
#'
#' @param Gradient an object returned by the function \code{\link{constructGradient}}
#' @param pred an object returned by applying the function \code{\link{predict}} to \code{Gradient}
#' @param measure whether plotted is species richness ("S"), an individual species ("Y") or community-weighted mean trait ("T")
#' @param index which species or trait is plotted
#' @param prob quantiles of the credibility interval plotted
#' @param showData whether raw data are plotted as well
#' @param jigger the amount by which the raw data are to be jiggered in x-direction (for factors) or y-direction (for continuous covariates)
#'
#' @return
#'
#' @details
#'
#' For \code{measure}="Y", \code{index} selects which species to plot from \code{hM$spNames}.
#' For \code{measure}="T", \code{index} selects which trait to plot from \code{hM$trNames}.
#' Whith \code{measure}="S" plotted is the row sum of \code{pred},
#' and thus the interpretation of "species richness" holds only for probit models.
#' For Poisson models "S" shows the total count,
#' whereas for normal models it shows the summed response.
#' For \code{measure}="T",
#' in probit model the weighting is over species occurrences,
#' whereas in count models it is over individuals.
#' In normal models, the weights are exp-transformed predictors to avoid negative weights
#'
#'
#' @seealso
#'
#' \code{\link{constructGradient}}, \code{\link{predict}}
#'
#' @examples
#'
#' Gradient = constructGradient(hM=m, focalVariable="x1")
#' predY = predict(m, Gradient=Gradient)
#' plotGradient(m, Gradient, pred=predY, measure="S")
#' plotGradient(m, Gradient, pred=predY, measure="Y", index = 2, showData = TRUE, jigger = 0.05)
#'
#' @export

plotGradient=function (hM, Gradient, predY, measure, index = 1, prob = c(0.025, 0.5, 0.975), showData = FALSE, jigger = 0,...){

   if (measure == "S") {
      predS = lapply(predY, rowSums)
      qpred = apply(abind(predS, along = 2), c(1), quantile,
                    prob = prob, na.rm=TRUE)
      ylabel = "Summed response"
      if (all(hM$distr[,1]==2)){
         ylabel = "Species richness"
      }
      if (all(hM$distr[,1]==3)){
         ylabel = "Total count"
      }
   }
   if (measure == "Y") {
      qpred = apply(abind(predY, along = 3), c(1, 2), quantile,
                    prob = prob, na.rm=TRUE)
      qpred = qpred[, , index]
      ylabel = hM$spNames[[index]]
   }
   if (measure == "T") {
      if (all(hM$distr[,1]==1)){
         predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)),
                                                                         hM$nt), ncol = hM$nt))
      } else {
         predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a),
                                                                    hM$nt), ncol = hM$nt))
      }
      qpred = apply(abind(predT, along = 3), c(1, 2), quantile,
                    prob = prob, na.rm = TRUE)
      qpred = qpred[, , index]
      ylabel = hM$trNames[[index]]
   }
   xlabel = colnames(Gradient$XDataNew)[[1]]
   xx = Gradient$XDataNew[, 1]
   lo = qpred[1, ]
   hi = qpred[3, ]
   me = qpred[2, ]

   lo1 = min(lo)
   hi1 = max(hi)

   if(showData){
      XDatacol = which(colnames(Gradient$XDataNew)[[1]]==colnames(hM$XData))
      if (measure == "S") {
         pY = rowSums(hM$Y)
      }
      if (measure == "Y") {
         pY = hM$Y[,index]
      }
      if (measure == "T") {
         if (all(hM$distr[,1]==1)){
            tmp = (exp(hM$Y) %*% hM$Tr)/matrix(rep(rowSums(exp(hM$Y)),hM$nt), ncol = hM$nt)

         } else {
            tmp = (hM$Y %*% hM$Tr)/matrix(rep(rowSums(hM$Y),hM$nt), ncol = hM$nt)
         }
         pY = tmp[,index]
      }
      pX = hM$XData[,XDatacol]

      lo1 = min(lo1,min(pY,na.rm = TRUE))
      hi1 = max(hi1,max(pY,na.rm = TRUE))
   }

   if (is.factor(xx)) {
      toPlot = data.frame(xx, me, lo, hi)
      plot.new()
      pl = ggplot(toPlot, aes_string(x = xx, y = me), fill = supp) + geom_bar(position = position_dodge(),
                                                                              stat = "identity") + xlab(xlabel) + ylab(ylabel) +
         geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2,
                       position = position_dodge(0.9))
      if(showData){
         if (jigger>0){
            pX = as.numeric(pX)
            pX = pX + runif(n =length(pY),min = -jigger, max = jigger)
         }
         dataToPlot = data.frame(pX = pX, pY = pY)
         pl = pl + geom_point(data = dataToPlot, aes_string(x = pX, y = pY))
      }
      plot(pl)
   } else {
      plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", xlab = xlabel, ylab = ylabel, ...)
      if(showData){
         if (jigger>0){
            de=(hi1-lo1)*jigger
            pY = lo1 + de + (hi1-lo1-2*de)*(pY-lo1)/(hi1-lo1) + runif(n =length(pY),min = -jigger, max = jigger)
         }
         dataToPlot = cbind(pX,pY)
         points(dataToPlot, pch = 16, col = "lightgrey")
      }
      polygon(c(xx, rev(xx)), c(qpred[1, ], rev(qpred[3, ])),
              col = rgb(0,0,1,alpha=.5), border = FALSE)
      lines(xx, qpred[2, ], lwd = 2)

   }
}
