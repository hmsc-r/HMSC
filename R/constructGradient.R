#' @title hM$sc$constructGradient
#'
#' @description Constructs an environmental gradient over one of the variables included in \code{XData}
#' @param focalVariable
#' @param non.focalVariable
#' @param ngrid
#'
#'
#' @return a named list with members \code{XDataNew}, \code{studyDesignNew} and \code{rLNew}
#'
#' @details
#'
#' non.focalVariables is a list, of which each element is of form variable=list(type,value).
#' where variable is one of the non-focal variables, and the value is needed only if type = 3
#' type = 1 sets the values of the non-focal variable
#' to the most likely value (defined as expected value for covarites, mode for factors)
#' type = 2 sets the values of the non-focal variable to most likely value, given the value of focal variable,
#' based on linear relationship
#' type = 3 fixes to the value given.
#' if a non.focalVariable is not listed, type=2 is used as default
#' note that if focal variable is continuous, selecting type 2 for a non-focal categorical variable can cause abrupt changes in response
#'
#' @seealso
#'
#' \code{\link{plotGradient}}
#'
#' @examples
#'
#' Gradient = constructGradient(hM=m, focalVariable="x2", non.focalVariables=list(x1=list(3,1),x3=list(1)))
#' predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew, ranLevels=Gradient$rLNew)
#' plotGradient(m, Gradient, pred=predY, measure="S")

#' @export

constructGradient = function(hM, focalVariable, non.focalVariables=list(), ngrid=20){

   non.focalNames = names(non.focalVariables)
   Mode <- function(x, na.rm=FALSE) {
      if(na.rm){
         x = x[!is.na(x)]
      }
      ux <- unique(x)
      return(ux[which.max(tabulate(match(x, ux)))])
   }

   vars = all.vars(hM$XFormula)
   nvars = length(vars)
   factors = rep(FALSE,nvars)
   focal = NA
   non.focals = NULL
   types = NULL
   vals = list()
   for (i in seq_len(nvars)){
      if (vars[i]==focalVariable){
         focal = i
      } else {
         non.focals = c(non.focals,i)
         found = FALSE
         for (j in seq_len(length(non.focalVariables))){
            if (vars[i]==non.focalNames[[j]]){
               found = TRUE
               type = as.numeric(non.focalVariables[[j]][[1]])
               types = c(types,type)
               if (type==3) {
                  vals[[length(vals)+1]] = non.focalVariables[[j]][[2]]
               } else{
                  vals[[length(vals)+1]] = NA
               }
            }
         }
         if (!found) {
            types = c(types,2)
            vals[[length(vals)+1]] = NA
         }
      }
      switch(class(hM$XData),
             matrix = {
                if (is.factor(hM$XData[,vars[i]])){
                   factors[i] = TRUE
                }
             },
             data.frame = {
                if (is.factor(hM$XData[,vars[i]])){
                   factors[i] = TRUE
                }
             },           list = {
                if (is.factor(hM$XData[[1]][,vars[i]])){
                   factors[i] = TRUE
                }
             }
      )
   }


   switch(class(hM$XData),
          matrix = {
             nz = 1
          },
          data.frame = {
             nz = 1
          },         list = {
             nz = hM$ns
             XDataNewList = list()
          }
   )

   for(case in 1:nz){
      switch(class(hM$XData),
             matrix = {
                XData = hM$XData
             },
             data.frame = {
                XData = hM$XData
             },           list = {
                XData = hM$XData[[case]]
             }
      )
      f.focal = factors[focal]
      v.focal = XData[,vars[focal]]
      if (f.focal){
         xx = levels(v.focal)
         ngrid = length(xx)
      } else {
         mi = min(v.focal)
         ma = max(v.focal)
         xx = seq(mi, ma, length.out = ngrid)
      }
      XDataNew = data.frame(xx)
      colnames(XDataNew) = vars[focal]

      for (i in seq_len(length(non.focals))){
         non.focal = non.focals[i]
         type = types[i]
         val = vals[[i]]
         f.non.focal = factors[non.focal]
         v.non.focal = XData[,vars[non.focal]]
         if (f.non.focal){
            if (type==1){
               XDataNew[,vars[non.focal]] = Mode(v.non.focal)
            }
            if (type==2){
               mymnm = multinom(v.non.focal~v.focal)
               yy=predict(mymnm, newdata=data.frame(v.focal=xx))
               XDataNew[,vars[non.focal]] = yy
            }
            if (type==3){
               XDataNew[,vars[non.focal]] = rep(val,ngrid)
            }
         }
         if (!f.non.focal){
            v.non.focal = XData[,vars[non.focal]]
            if (type==1){
               XDataNew[,vars[non.focal]] = mean(v.non.focal)
            }
            if (type==2){
               mylm = lm(v.non.focal~v.focal)
               yy=predict(mylm, newdata=data.frame(v.focal=xx))
               XDataNew[,vars[non.focal]] = yy
            }
            if (type==3){
               XDataNew[,vars[non.focal]] = rep(val,ngrid)
            }
         }
      }
      if(class(hM$XData)=="list"){XDataNewList[[case]]=XDataNew}
   }
   if(class(hM$XData)=="list"){XDataNew = XDataNewList}

   dfPiNew = matrix(NA,ngrid,hM$nr)
   colnames(dfPiNew) = hM$rLNames
   for (r in seq_len(hM$nr)){
      dfPiNew[,r] = sprintf('new_unit',1:(ngrid))
   }
   dfPiNew = as.data.frame(dfPiNew)

   rLNew = vector("list", hM$nr)
   names(rLNew) = hM$rLNames
   for (r in seq_len(hM$nr)){
      rL1 = hM$rL[[r]]
      units = rL1$pi
      units1 = c(units,"new_unit")
      xydata = rL1$s
      if (!is.null(xydata)){
         nxy = dim(xydata)[1]
         xydata1 = matrix(NA,nrow = nxy+1, ncol = rL1$sDim)
         xydata1[1:nxy,] = xydata
         xydata1[nxy+1,] = colMeans(xydata)
         rownames(xydata1) = units1
         colnames(xydata1) = colnames(xydata)
         rL1$s = xydata1
      }
      rL1$pi = units1
      rL1$N = rL1$N+1
      rLNew[[r]] = rL1
   }

   Gradient = list(XDataNew=XDataNew, studyDesignNew=dfPiNew, rLNew=rLNew)
   return(Gradient)
}
