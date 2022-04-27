#' @export

print.Hmsc = function(x, ...)
{
    cat(sprintf(
        "Hmsc object with %d sampling units, %d species, %d covariates, %d traits and %d random levels\n",
        x$ny, x$ns, x$nc, x$nt, x$nr))
    ## Information about sampleMcmc if done on the same object
    if (!is.null(x$samples) && x$samples > 0) {
        cat("Posterior MCMC sampling with ")
        if ((nchain <- length(x$postList)) > 1 )
            cat(nchain, "chains each with ")
        cat(x$samples, "samples, thin", x$thin, "and transient", x$transient,
            "\n")
        ## info on failed updaters (if any)
        fails <- sapply(x$postList, attr, which = "failedUpdates")
        if (is.matrix(fails) && any(fails > 0)) {
            fails <- t(fails)
            rownames(fails) <- paste0("Chain ",
                                      as.character(seq_len(nrow(fails))), ":")
            fails <- fails[rowSums(fails) > 0, colSums(fails) > 0, drop = FALSE]
            cat("Updater failures in following chains with ",
                x$samples * x$thin + x$transient, " attemps in each chain:\n")
            print(fails)
        }
    } else {
        cat("No posterior samples\n")
    }
}
