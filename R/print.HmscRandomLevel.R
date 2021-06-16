#' @export

print.HmscRandomLevel = function(x, ...){
   cat(sprintf("Hmsc random level object with %d units. Number of covariates is %d.\n",
               x$N, x$xDim))
}
