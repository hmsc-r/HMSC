#' @export

print.Hmsc = function(x, ...)
{
    cat(sprintf(
        "Hmsc object with %d sampling units, %d species, %d covariates, %d traits and %d random levels\n",
        x$ny, x$ns, x$nc, x$nt, x$nr))
    ## Information about sampleMcmc if done on the same object
    if (!is.null(x$sample) && x$sample > 0) {
        cat("Posterior MCMC sampling with ")
        if ((nchain <- length(x$postList)) > 1 )
            cat(nchain, "chains each with ")
        cat(x$samples, "samples, thin", x$thin, "and transient", x$transient,
            "\n")
    }
}
