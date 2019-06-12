#' @export

print.Hmsc = function(x, ...){
   cat(sprintf("Hmsc object with %d sampling units, %d species, %d covariates, %d traits and %d random levels\n", x$ny, x$ns, x$nc, x$nt, x$nr))
}
