#' @title plotGradient
#'
#' @description Plots an environmental gradient over one of the variables included in \code{XData}
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param Gradient an object returned by \code{\link{constructGradient}}
#' @param predY an object returned by applying the function \code{\link{predict}} to \code{Gradient}
#' @param measure whether to plot species richness ("S"), an individual species ("Y") or community-weighted
#' mean trait values ("T")
#' @param index which species or trait to plot
#' @param q quantiles of the credibility interval plotted
#' @param xlabel label for x-axis
#' @param ylabel label for y-axis
#' @param cicol colour with which the credibility interval is plotted
#' @param pointcol colour with which the data points are plotted
#' @param pointsize size in which the data points are plotted
#' @param showData whether raw data are plotted as well
#' @param jigger the amount by which the raw data are to be jiggered in x-direction (for factors) or
#' y-direction (for continuous covariates)
#'
#' @param yshow scale y-axis so that these values are also
#'     visible. This can used to scale y-axis so that it includes 0
#'     and the expected maximum values.
#'
#' @param ... additional arguments for plot
#'
#' @return For the case of a continuous covariate, returns the posterior probability that the plotted
#' variable is greater for the last sampling unit of the gradient than for the first sampling unit of
#' the gradient. For the case of a factor, returns the plot object.
#'
#' @details
#'
#' For \code{measure}="Y", \code{index} selects which species to plot from \code{hM$spNames}.
#' For \code{measure}="T", \code{index} selects which trait to plot from \code{hM$trNames}.
#' With \code{measure}="S" the row sum of \code{pred} is plotted,
#' and thus the interpretation of "species richness" holds only for probit models.
#' For Poisson models "S" shows the total count,
#' whereas for normal models it shows the summed response.
#' For \code{measure}="T",
#' in probit model the weighting is over species occurrences,
#' whereas in count models it is over individuals.
#' In normal models, the weights are exp-transformed predictions to avoid negative weights
#'
#'
#' @seealso
#' \code{\link{constructGradient}}, \code{\link{predict}}
#'
#' @examples
#' # Plot response of species 2 over the gradient of environmental variable x1
#' Gradient = constructGradient(TD$m, focalVariable="x1")
#' predY = predict(TD$m, Gradient=Gradient)
#' plotGradient(TD$m, Gradient, pred=predY, measure="Y", index = 2, showData = TRUE, jigger = 0.05)
#' # Plot modelled species richness over the gradient of environmental variable x1
#' Gradient = constructGradient(TD$m, focalVariable="x1")
#' predY = predict(TD$m, Gradient=Gradient)
#' plotGradient(TD$m, Gradient, pred=predY, measure="S")
#'
#' @importFrom stats quantile runif
#' @importFrom graphics plot axis points polygon lines
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot aes_string geom_bar position_dodge xlab ylab geom_errorbar
#'   aes geom_point
#' @importFrom abind abind
#'
#' @export

plotGradient =
    function (hM, Gradient, predY, measure, xlabel = NULL, ylabel = NULL,
              index = 1, q = c(0.025, 0.5, 0.975), cicol = rgb(0,0,1,alpha=.5),
              pointcol = "lightgrey", pointsize = 1, showData = FALSE,
              jigger = 0, yshow = NA,  ...)
{

   Pr = NA

   if(is.null(xlabel)){
      switch(class(hM$X)[1L],
             matrix = {
                xlabel = colnames(Gradient$XDataNew)[[1]]
             },
             list = {
                xlabel = colnames(Gradient$XDataNew[[1]])[[1]]
             }
      )
   }

   switch(class(hM$X)[1L],
          matrix = {
             xx = Gradient$XDataNew[, 1]
          },
          list = {
             if (measure == "Y") {
                xx = Gradient$XDataNew[[index]][, 1]
             } else {
                xx = Gradient$XDataNew[[1]][, 1]
             }
          }
   )

   ngrid = length(xx)

   if (measure == "S") {
      predS = abind(lapply(predY, rowSums),along=2)
      Pr = mean(predS[ngrid,]>predS[1,])
      qpred = apply(predS, c(1), quantile,
                    probs = q, na.rm=TRUE)
      if(is.null(ylabel)){
         ylabel = "Summed response"
         if (all(hM$distr[,1]==2)){
            ylabel = "Species richness"
         }
         if (all(hM$distr[,1]==3)){
            ylabel = "Total count"
         }
      }
   }
   if (measure == "Y") {
      tmp = abind(predY, along = 3)
      Pr = mean(tmp[ngrid,index,]>tmp[1,index,])
      qpred = apply(tmp, c(1, 2), quantile,
                    probs = q, na.rm=TRUE)
      qpred = qpred[, , index]
      if(is.null(ylabel)){
         ylabel = hM$spNames[[index]]
      }
   }
   if (measure == "T") {
      if (all(hM$distr[,1]==1)){
         predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)),
                                                                         hM$nt), ncol = hM$nt))
      } else {
         predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a),
                                                                    hM$nt), ncol = hM$nt))
      }
      predT = abind(predT, along = 3)
      Pr = mean(predT[ngrid,index,]>predT[1,index,])
      qpred = apply(predT, c(1, 2), quantile,
                    probs = q, na.rm = TRUE)
      qpred = qpred[, , index]
      if(is.null(ylabel)){
         ylabel = hM$trNames[[index]]
      }
   }

   lo = qpred[1, ]
   hi = qpred[3, ]
   me = qpred[2, ]

   lo1 = min(lo, yshow, na.rm = TRUE)
   hi1 = max(hi, yshow, na.rm = TRUE)

   if(showData){
      switch(class(hM$X)[1L],
             matrix = {
                XDatacol = which(colnames(Gradient$XDataNew)[[1]]==colnames(hM$XData))
             },
             list = {
                XDatacol = which(colnames(Gradient$XDataNew[[1]])[[1]]==colnames(hM$XData[[1]]))
             }
      )
      if (measure == "S") {
         pY = rowSums(hM$Y, na.rm = TRUE)
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
      if (!is.numeric(pY))
         pY <- as.numeric(pY)

      switch(class(hM$X)[1L],
             matrix = {
                pX = hM$XData[,XDatacol]
             },
             list = {
                pX = hM$XData[[1]][,XDatacol]
             }
      )

      hi1 <- max(hi1, max(pY, na.rm = TRUE))
      lo1 <- min(lo1, min(pY, na.rm = TRUE))
   }

   if (is.factor(xx)) {
      toPlot = data.frame(xx, me, lo, hi, stringsAsFactors = TRUE)
#      plot.new()
      pl = ggplot(toPlot, aes_string(x = xx, y = me)) + geom_bar(position = position_dodge(),
                                                                              stat = "identity") + xlab(xlabel) + ylab(ylabel) +
         geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2,
                       position = position_dodge(0.9))
      if(showData){
         if (jigger>0){
            pX = as.numeric(pX)
            pX = pX + runif(n =length(pY),min = -jigger, max = jigger)
         }
         dataToPlot = data.frame(pX = pX, pY = pY, stringsAsFactors = TRUE)
         pl = pl + geom_point(data = dataToPlot, aes_string(x = pX, y = pY), size = pointsize)
      }
      #plot(pl)
   } else {
      if (inherits(hM$X,"list") && !measure=="Y") {
         plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", xaxt = "n", xlab = xlabel, ylab = ylabel, ...)
         axis(1,c(min(xx),(min(xx)+max(xx))/2,max(xx)),c("min","mean","max"))
      } else {
         plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", xlab = xlabel, ylab = ylabel, ...)
      }
      if(showData){
         if (jigger>0){
            de=(hi1-lo1)*jigger
            pY = lo1 + de + (hi1-lo1-2*de)*(pY-lo1)/(hi1-lo1) + runif(n =length(pY),min = -jigger, max = jigger)
         }
         dataToPlot = cbind(pX,pY)
         points(dataToPlot, pch = 16, col = pointcol)
      }
      polygon(c(xx, rev(xx)), c(qpred[1, ], rev(qpred[3, ])),
              col =  cicol, border = FALSE)
      lines(xx, qpred[2, ], lwd = 2)

   }
   return(if(is.factor(xx)) pl else Pr)
}
