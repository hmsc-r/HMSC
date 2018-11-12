#' @title Hmsc$plotGradient
#'
#' @description Plots ...
#'
#' @param Gradient
#' @param pred
#' @param measure
#' @param index
#' @param prob
#'
#'
#' @return
#'
#'
#' @seealso
#'
#' 
#' @examples
#'

plotGradient = function(Gradient, pred, measure, index=1, prob=c(0.025,0.5,0.975)){

  if (measure == "S"){
    predS = lapply(predY, rowSums)
    qpred = apply(abind(predS,along=2),c(1),quantile,prob=prob)
    ylabel = "Species richness"
  }
  if (measure == "Y"){
    qpred = apply(abind(predY,along=3),c(1,2),quantile,prob=prob)
    qpred = qpred[,,index]
    ylabel = self$spNames[[index]]
  }
  if (measure == "T"){
    predT = lapply(predY, function(a) (a%*%self$Tr)/matrix(rep(rowSums(a),self$nt),ncol=self$nt))
    qpred = apply(abind(predT,along=3),c(1,2),quantile,prob=prob)
    qpred = qpred[,,index]
    ylabel = self$trNames[[index]]
  }

  xlabel = colnames(Gradient$XDataNew)[[1]]
  xx = Gradient$XDataNew[,1]
  if (is.factor(xx)){
    lo = qpred[1,]
    hi = qpred[3,]
    me = qpred[2,]
    toPlot = data.frame(xx,me,lo,hi)
    ggplot(toPlot, aes(x=xx, y=me),fill=supp) +
      geom_bar(position=position_dodge(), stat="identity") +
      xlab(xlabel) + ylab(ylabel) +
      geom_errorbar(aes(ymin=lo, ymax=hi),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))
  } else
  {
    plot(xx, qpred[2,], ylim = c(min(qpred[1,]),max(qpred[3,])), type = "l", xlab = xlabel, ylab = ylabel)
    polygon(c(xx,rev(xx)),c(qpred[1,],rev(qpred[3,])),col = "grey75", border = FALSE)
    lines(xx,qpred[2,], lwd = 2)
    lines(xx, qpred[1,], col="red",lty=2)
    lines(xx, qpred[3,], col="red",lty=2)
  }
}

Hmsc$set("public", "plotGradient", plotGradient, overwrite=TRUE)

