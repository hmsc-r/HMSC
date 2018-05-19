#' @title Hmsc$constructGradient
#'
#' @description Computes ...
#' @param focalVariable
#' @param non.focalVariable
#' @param ngrid
#'
#' @examples
#'

constructGradient = function(focalVariable, non.focalVariables=list(), ngrid=20){

  non.focalNames = names(non.focalVariables)
  Mode <- function(x, na.rm = FALSE) {
    if(na.rm){
      x = x[!is.na(x)]
    }
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }

  vars = all.vars(self$XFormula)
  nvars = length(vars)
  factors = rep(FALSE,nvars)
  focal = NA
  non.focals = NULL
  types = NULL
  vals = list()
  for (i in 1:nvars){
    if (vars[i]==focalVariable){
      focal = i
    } else {
      non.focals = c(non.focals,i)
      found = FALSE
      if (length(non.focalVariables)>0){
        for (j in 1:length(non.focalVariables)){
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
      }
      if (!found) {
        types = c(types,2)
        vals[[length(vals)+1]] = NA
      }
    }
    if (is.factor(self$XData[,vars[i]])){
      factors[i] = TRUE
    }
  }


  f.focal = factors[focal]
  v.focal = self$XData[,vars[focal]]
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

  for (i in 1:length(non.focals)){
    non.focal = non.focals[i]
    type = types[i]
    val = vals[[i]]
    f.non.focal = factors[non.focal]
    v.non.focal = self$XData[,vars[non.focal]]
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
      v.non.focal = self$XData[,vars[non.focal]]
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


  dfPiNew = matrix(NA,ngrid,self$nr)
  for (r in 1:self$nr){
    dfPiNew[,r] = sprintf('new_unit',1:(ngrid))
  }
  dfPiNew = as.data.frame(dfPiNew)

  rLNew = vector("list", self$nr)
  for (r in 1:self$nr){
    tmp = self$rL[[r]]
    rL1=tmp$clone()
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

  Gradient = list(XDataNew = XDataNew, dfPiNew = dfPiNew, rLNew = rLNew)

  return(Gradient)
}

Hmsc$set("public", "constructGradient", constructGradient, overwrite=TRUE)
