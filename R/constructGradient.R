#' @title constructGradient
#'
#' @description Constructs an environmental gradient over one of the variables included in \code{XData}
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param focalVariable focal variable over which the gradient is constructed
#' @param non.focalVariables list giving assumptions on how non-focal variables co-vary with the focal variable or a single number given the default type for all non-focal variables
#' @param ngrid number of points along the gradient (for continuous focal variables)
#'
#' @return a named list with slots \code{XDataNew}, \code{studyDesignNew} and \code{rLNew}
#'
#' @details
#' In basic form, \code{non.focalVariables} is a list, where each element is on the form variable=list(type,value),
#' where \code{variable} is one of the non-focal variables, and the \code{value} is needed only if \code{type = 3}. Alternatives
#' \code{type = 1} sets the values of the non-focal variable
#' to the most likely value (defined as expected value for covariates, mode for factors),
#' \code{type = 2} sets the values of the non-focal variable to most likely value, given the value of focal variable,
#' based on a linear relationship, and
#' \code{type = 3} fixes to the value given.
#' As a shortcut, a single number \code{1} or \code{2} can be given as a type
#' used for all non-focal variables.
#' If a \code{non.focalVariable} is not listed, \code{type=2} is used as default.
#' Note that if the focal variable is continuous, selecting type 2 for a non-focal categorical variable can cause abrupt changes in response.
#'
#' @seealso
#' \code{\link{plotGradient}}, \code{\link{predict}}.
#'
#' @examples
#' # Construct gradient for environmental covariate called 'x1'.
#' Gradient = constructGradient(TD$m, focalVariable="x1")
#'
#' # Construct gradient for environmental covariate called 'x1'
#' # while setting the other covariate to its most likely values
#' Gradient = constructGradient(TD$m, focalVariable="x1",non.focalVariables=list(x2=list(1)))
#'
#' @importFrom stats lm predict
#' @importFrom methods is
#' @importFrom sp coordinates `coordinates<-` proj4string `proj4string<-`
#' @importFrom nnet multinom
#'
#' @export

constructGradient =
   function(hM, focalVariable, non.focalVariables=list(), ngrid=20)
{
   ## default type 2 unless a single number is given as a non-focal variable
   if (is.numeric(non.focalVariables) && length(non.focalVariables) == 1) {
      defType <- non.focalVariables
      if (!(defType %in% 1:2)) {
         warning("only type 1 and 2 are allowed as a single number: setting to 2")
         defType <- 2
      }
      non.focalVariables <- NULL
   } else {
      defType <- 2
   }

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
            types = c(types, defType)
            vals[[length(vals)+1]] = NA
         }
      }
      switch(class(hM$XData)[1L],
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


   switch(class(hM$XData)[1L],
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
      switch(class(hM$XData)[1L],
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
      XDataNew = data.frame(xx, stringsAsFactors = TRUE)
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
               yy=predict(mymnm,
                          newdata=data.frame(v.focal=xx,
                                             stringsAsFactors = TRUE))
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
               yy=predict(mylm,
                          newdata=data.frame(v.focal=xx,
                                             stringsAsFactors = TRUE))
               XDataNew[,vars[non.focal]] = yy
            }
            if (type==3){
               XDataNew[,vars[non.focal]] = rep(val,ngrid)
            }
         }
      }
      if(inherits(hM$XData, "list")) {XDataNewList[[case]]=XDataNew}
   }
   if(inherits(hM$XData, "list")){XDataNew = XDataNewList}

   dfPiNew = matrix("new_unit", ngrid, hM$nr)
   colnames(dfPiNew) = hM$rLNames
   dfPiNew = as.data.frame(dfPiNew, stringsAsFactors = TRUE)

   rLNew = vector("list", hM$nr)
   names(rLNew) = hM$rLNames
   for (r in seq_len(hM$nr)){
      rL1 = hM$rL[[r]]
      xydata = rL1$s
      if (!is.null(xydata)) {
         if (is(xydata, "Spatial")) {
            centre <- as.data.frame(t(colMeans(coordinates(xydata))))
            rownames(centre) <- "new_unit"
            coordinates(centre) <- colnames(centre)
            proj4string(centre) <- proj4string(xydata)
            xydata <- rbind(xydata, centre)
         } else {
            xydata = rbind(xydata, "new_unit" = colMeans(xydata))
         }
      } # end !is.null(xydata)
      rL1$s = xydata
      distMat = rL1$distMat
      if (!is.null(distMat)){
         units1 = c(rownames(distMat), "new_unit")
         rm=rowMeans(distMat)
         focals = order(rm)[1:2]
         newdist = colMeans(distMat[focals,])
         distMat1=cbind(distMat,newdist)
         distMat1=rbind(distMat1,c(newdist,0))
         rownames(distMat1) = units1
         colnames(distMat1) = units1
         rL1$distMat = distMat1
      }
      rL1$pi = c(rL1$pi, "new_unit")
      rL1$N = rL1$N+1
      rLNew[[r]] = rL1
   }

   Gradient = list(XDataNew=XDataNew, studyDesignNew=dfPiNew, rLNew=rLNew)
   return(Gradient)
}
