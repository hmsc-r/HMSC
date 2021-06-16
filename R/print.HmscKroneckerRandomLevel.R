#' @export

print.HmscKroneckerRandomLevel = function(x, ...){
   cat(sprintf("Kronecker Hmsc random level object made from %d Hmsc random levels with %d unit combinations",
               length(x$rLList), x$N))
}
