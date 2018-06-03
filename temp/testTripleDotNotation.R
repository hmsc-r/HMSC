plotVariancePartitioning = function(...){
   mainTitle = "barplot"
   ng = 1
   a = list(...)
   par(...)
   barplot(1:10, main=mainTitle, xlab="Species", ylab="Variance proportion", col=heat.colors(ng,alpha=1), ...)
}

plotVariancePartitioning(main="aaa")
