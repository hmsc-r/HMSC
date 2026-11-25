check_fitted_postprocessing <- function(fitted_model, nSamples, nChains) {
   coda_post <- convertToCodaObject(fitted_model)

   expect_true(coda::is.mcmc.list(coda_post$Beta))
   expect_equal(length(coda_post$Beta), nChains)
   expect_equal(nrow(coda_post$Beta[[1]]), nSamples)
   expect_equal(ncol(coda_post$Beta[[1]]), fitted_model$nc * fitted_model$ns)

   beta_estimate <- getPostEstimate(fitted_model, parName = "Beta", q = c(0.025, 0.975))
   expect_equal(dim(beta_estimate$q), c(2, fitted_model$nc, fitted_model$ns))

   sigma_matrix <- do.call(rbind, lapply(fitted_model$postList, function(chain_samples) {
      do.call(rbind, lapply(chain_samples, function(posterior_sample) {
         posterior_sample$sigma
      }))
   }))
   sigma_estimate <- getPostEstimate(fitted_model, parName = "sigma", q = c(0.025, 0.975))
   expect_equal(dim(sigma_estimate$q), c(2, fitted_model$ns))
   expect_equal(
      unname(sigma_estimate$q),
      unname(apply(sigma_matrix, 2, stats::quantile, c(0.025, 0.975)))
   )

   if (fitted_model$phyloFlag) {
      expect_true(coda::is.mcmc.list(coda_post$rho))
      expect_equal(length(coda_post$rho), nChains)
      rho_matrix <- do.call(rbind, coda_post$rho)
      expected_rho <- do.call(rbind, lapply(fitted_model$postList, function(chain_samples) {
         do.call(rbind, lapply(chain_samples, function(posterior_sample) {
            fitted_model$rhopw[posterior_sample$rhoInd, 1]
         }))
      }))
      expect_equal(unname(rho_matrix), unname(expected_rho))

      rho_estimate <- getPostEstimate(fitted_model, parName = "rho", q = c(0.025, 0.975))
      expect_equal(dim(rho_estimate$q), c(2, fitted_model$rhoLen))
      expect_equal(
         unname(rho_estimate$q),
         unname(apply(rho_matrix, 2, stats::quantile, c(0.025, 0.975)))
      )
      expect_null(convertToCodaObject(fitted_model, rho = FALSE)$rho)
   } else {
      expect_null(coda_post$rho)
      expect_warning(
         coda_forced_rho <- convertToCodaObject(fitted_model, rho = TRUE),
         "without phylogeny"
      )
      expect_null(coda_forced_rho$rho)
   }

   for (chain_idx in seq_along(fitted_model$postList)) {
      for (sample_idx in seq_along(fitted_model$postList[[chain_idx]])) {
         posterior_sample <- fitted_model$postList[[chain_idx]][[sample_idx]]
         expect_equal(length(posterior_sample$alpha), fitted_model$nr)
         expect_equal(length(posterior_sample$alphaInd), fitted_model$nr)
         for (level_idx in seq_len(fitted_model$nr)) {
            expect_null(dim(posterior_sample$alphaInd[[level_idx]]))
            expect_equal(
               length(posterior_sample$alphaInd[[level_idx]]),
               ncol(posterior_sample$Eta[[level_idx]])
            )
            if (fitted_model$rL[[level_idx]]$sDim > 0) {
               expect_null(dim(posterior_sample$alpha[[level_idx]]))
            } else {
               expect_null(posterior_sample$alpha[[level_idx]])
            }
         }
      }
   }

   if (fitted_model$nr == 0) {
      expect_null(coda_post$alpha)
      expect_null(coda_post$Eta)
      expect_warning(
         coda_forced_alpha <- convertToCodaObject(fitted_model, alpha = TRUE),
         "without spatial random levels"
      )
      expect_null(coda_forced_alpha$alpha)
   } else {
      has_spatial <- any(vapply(fitted_model$rL, function(rL) rL$sDim > 0, logical(1)))
      if (has_spatial) {
         expect_equal(length(coda_post$alpha), fitted_model$nr)
      } else {
         expect_null(coda_post$alpha)
         coda_nonspatial_defaults <- convertToCodaObject(
            fitted_model, Eta = FALSE, Lambda = FALSE, Omega = FALSE,
            Psi = FALSE, Delta = FALSE
         )
         expect_null(coda_nonspatial_defaults$alpha)
         expect_warning(
            coda_forced_alpha <- convertToCodaObject(fitted_model, alpha = TRUE),
            "without spatial random levels"
         )
         expect_null(coda_forced_alpha$alpha)
      }

      for (level_idx in seq_len(fitted_model$nr)) {
         if (fitted_model$rL[[level_idx]]$sDim > 0) {
            expect_true(coda::is.mcmc.list(coda_post$alpha[[level_idx]]))
            expect_equal(length(coda_post$alpha[[level_idx]]), nChains)
            expect_equal(nrow(coda_post$alpha[[level_idx]][[1]]), nSamples)

            alpha_matrix <- do.call(rbind, coda_post$alpha[[level_idx]])
            expected_alpha <- do.call(rbind, lapply(fitted_model$postList, function(chain_samples) {
               do.call(rbind, lapply(chain_samples, function(posterior_sample) {
                  fitted_model$rL[[level_idx]]$alphapw[posterior_sample$alphaInd[[level_idx]], 1]
               }))
            }))
            expect_equal(unname(alpha_matrix), unname(expected_alpha))
            expect_true(all(alpha_matrix %in% fitted_model$rL[[level_idx]]$alphapw[, 1]))

            alpha_estimate <- getPostEstimate(
               fitted_model, parName = "alpha", r = level_idx, q = c(0.025, 0.5, 0.975)
            )
            expect_equal(unname(alpha_estimate$mean), unname(colMeans(alpha_matrix)))
            expect_equal(dim(alpha_estimate$q), c(3, ncol(alpha_matrix)))
            expect_equal(
               unname(alpha_estimate$q),
               unname(apply(alpha_matrix, 2, stats::quantile, c(0.025, 0.5, 0.975)))
            )

            alpha_median <- getPostEstimate(
               fitted_model, parName = "alpha", r = level_idx, q = 0.5
            )
            expect_equal(length(alpha_median$q), ncol(alpha_matrix))
            expect_equal(
               unname(alpha_median$q),
               unname(apply(alpha_matrix, 2, stats::quantile, 0.5))
            )
         } else {
            if (has_spatial) {
               expect_null(coda_post$alpha[[level_idx]])
            }
            expect_error(
               getPostEstimate(fitted_model, parName = "alpha", r = level_idx),
               "is not spatial and has no alpha parameter"
            )
         }
      }

      coda_no_alpha <- convertToCodaObject(fitted_model, alpha = FALSE)
      expect_null(coda_no_alpha$alpha)
      expect_true(coda::is.mcmc.list(coda_no_alpha$Beta))

      coda_lambda_only <- if (has_spatial) {
         convertToCodaObject(fitted_model, Eta = FALSE, Lambda = TRUE, alpha = TRUE)
      } else {
         convertToCodaObject(fitted_model, Eta = FALSE, Lambda = TRUE)
      }
      expect_equal(coda_lambda_only$alpha, coda_post$alpha)

      coda_fixed_only <- convertToCodaObject(
         fitted_model, Eta = FALSE, Lambda = FALSE, Omega = FALSE,
         alpha = FALSE, Psi = FALSE, Delta = FALSE
      )
      expect_true(coda::is.mcmc.list(coda_fixed_only$Beta))
      expect_null(coda_fixed_only$alpha)
   }

   variance_partitioning <- computeVariancePartitioning(fitted_model)
   expect_true(is.matrix(variance_partitioning$vals))
   expect_equal(ncol(variance_partitioning$vals), fitted_model$ns)

   predictions <- predict(fitted_model, expected = TRUE)
   expect_equal(length(predictions), nChains * nSamples)
   expect_identical(dim(predictions[[1]]), c(fitted_model$ny, fitted_model$ns))
}

# Specific mimic helper designed solely to plug the Hmsc-HPC gap in pure-R tests
# by reshaping sampleMcmc(..., engine = "HPC")$initParList into the posterior
# sample list layout expected by importPosteriorFromHPC().
# Note: Unless there is super solid evidence that this specific mimic piece is
# incorrect, there is no need to consult the whole Hmsc-HPC codebase. If
# specifically this mimic piece requires a fix, then the up-to-date Hmsc-HPC
# codebase shall be investigated.
convert_init_par_to_hpc_chains <- function(hpc_init_object, model, nSamples) {
   lapply(hpc_init_object$initParList, function(init_par) {
      lapply(seq_len(nSamples), function(sample_idx) {
         posterior_sample <- init_par
         posterior_sample$iV <- chol2inv(chol(init_par$V))
         posterior_sample$V <- NULL
         posterior_sample$Z <- NULL
         posterior_sample$sigma <- sqrt(init_par$sigma)
         if (model$phyloFlag) {
            posterior_sample$rhoInd <- rep(
               ((sample_idx - 1L) %% nrow(model$rhopw)) + 1L,
               model$rhoLen
            )
         }
         for (level_idx in seq_len(model$nr)) {
            if (model$rL[[level_idx]]$sDim > 0) {
               alpha_grid_size <- nrow(model$rL[[level_idx]]$alphapw)
               posterior_sample$alphaInd[[level_idx]] <-
                  ((seq_along(init_par$alphaInd[[level_idx]]) + sample_idx - 1L) %% alpha_grid_size) + 1L
            }
         }
         posterior_sample
      })
   })
}

test_that("R fitting and post-processing work across random-level structures", {
   set.seed(123)
   nSamples <- 5
   transient <- 2
   thin <- 1
   nChains <- 2

   knots <- constructKnots(TD$xycoords, knotDist = 0.5, minKnotDist = 0.3)
   rL_spatial_gpp <- HmscRandomLevel(sData = TD$xycoords, sMethod = "GPP", sKnot = knots)
   rL_spatial_nngp <- HmscRandomLevel(sData = TD$xycoords, sMethod = "NNGP", nNeighbours = 4)
   rL_plot_nonspatial <- HmscRandomLevel(units = levels(TD$studyDesign$plot))

   models <- list(
      no_ranLevels = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2
      ),
      single_nonspatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(sample = TD$rL2)
      ),
      single_spatial_full = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(plot = TD$rL1)
      ),
      single_spatial_gpp = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(plot = rL_spatial_gpp)
      ),
      single_spatial_nngp = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(plot = rL_spatial_nngp)
      ),
      spatial_then_nonspatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         phyloTree = TD$phy,
         studyDesign = TD$studyDesign, ranLevels = list(plot = TD$rL1, sample = TD$rL2)
      ),
      nonspatial_then_spatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(sample = TD$rL2, plot = TD$rL1)
      ),
      two_nonspatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(sample = TD$rL2, plot = rL_plot_nonspatial)
      )
   )

   for (model_name in names(models)) {
      fitted_model <- sampleMcmc(
         models[[model_name]], samples = nSamples, transient = transient,
         thin = thin, nChains = nChains, verbose = 0
      )
      check_fitted_postprocessing(fitted_model, nSamples = nSamples, nChains = nChains)
   }

   prior_fitted_model <- sampleMcmc(
      models$spatial_then_nonspatial, samples = nSamples,
      nChains = nChains, fromPrior = TRUE
   )
   check_fitted_postprocessing(prior_fitted_model, nSamples = nSamples, nChains = nChains)

   missing_spatial_alpha_model <- prior_fitted_model
   missing_spatial_alpha_model$postList[[1]] <- lapply(
      missing_spatial_alpha_model$postList[[1]],
      function(posterior_sample) {
         posterior_sample$alpha[1] <- list(NULL)
         posterior_sample
      }
   )
   expect_error(convertToCodaObject(missing_spatial_alpha_model))
})

test_that("HPC export and importPosteriorFromHPC round-trip preserves alpha and post-processing", {
   set.seed(456)
   nSamples <- 5
   transient <- 2
   thin <- 1
   nChains <- 2

   knots <- constructKnots(TD$xycoords, knotDist = 0.5, minKnotDist = 0.3)
   rL_spatial_gpp <- HmscRandomLevel(sData = TD$xycoords, sMethod = "GPP", sKnot = knots)

   models <- list(
      single_nonspatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(sample = TD$rL2)
      ),
      single_spatial_gpp = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         studyDesign = TD$studyDesign, ranLevels = list(plot = rL_spatial_gpp)
      ),
      spatial_then_nonspatial = Hmsc(
         Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
         TrData = TD$Tr, TrFormula = ~T1 + T2,
         phyloTree = TD$phy,
         studyDesign = TD$studyDesign, ranLevels = list(plot = TD$rL1, sample = TD$rL2)
      )
   )

   for (model_name in names(models)) {
      model <- models[[model_name]]
      hpc_init_object <- sampleMcmc(
         model, samples = nSamples, transient = transient,
         thin = thin, nChains = nChains, verbose = 0, engine = "HPC"
      )
      hpc_chains <- convert_init_par_to_hpc_chains(hpc_init_object, model, nSamples)
      hpc_fitted_model <- importPosteriorFromHPC(
         model, hpc_chains, nSamples = nSamples, thin = thin, transient = transient
      )
      check_fitted_postprocessing(hpc_fitted_model, nSamples = nSamples, nChains = nChains)

      coda_post <- convertToCodaObject(hpc_fitted_model)
      for (level_idx in seq_len(hpc_fitted_model$nr)) {
         if (hpc_fitted_model$rL[[level_idx]]$sDim > 0) {
            alpha_matrix <- do.call(rbind, coda_post$alpha[[level_idx]])
            expect_gt(max(alpha_matrix), 0)
            expect_gt(length(unique(as.vector(alpha_matrix))), 1L)
         }
      }
      if (hpc_fitted_model$phyloFlag) {
         rho_matrix <- do.call(rbind, coda_post$rho)
         expect_gt(max(rho_matrix), 0)
         expect_gt(length(unique(as.vector(rho_matrix))), 1L)
      }

      if (model_name %in% c("single_nonspatial", "spatial_then_nonspatial")) {
         truncated_chains <- lapply(hpc_chains, function(chain_samples) {
            lapply(chain_samples, function(posterior_sample) {
               if (model_name == "single_nonspatial") {
                  posterior_sample$alphaInd <- list()
               } else {
                  posterior_sample$alphaInd[[2]] <- NULL
               }
               posterior_sample
            })
         })
         truncated_fit <- importPosteriorFromHPC(
            model, truncated_chains, nSamples = nSamples, thin = thin, transient = transient
         )
         check_fitted_postprocessing(truncated_fit, nSamples = nSamples, nChains = nChains)
      }
   }
})

test_that("alignPosterior keeps alpha aligned with alphaInd across chains with differing nf", {
   set.seed(789)
   rL_spatial_custom <- setPriors(
      TD$rL1,
      alphapw = cbind(c(0.1, 0.5, 1.0), c(0.4, 0.3, 0.3)),
      nfMin = 2, nfMax = 4
   )
   rL_nonspatial_custom <- setPriors(TD$rL2, nfMin = 2, nfMax = 4)
   model <- Hmsc(
      Y = TD$Y, XData = TD$X, XFormula = ~x1 + x2, distr = "probit",
      TrData = TD$Tr, TrFormula = ~T1 + T2,
      studyDesign = TD$studyDesign,
      ranLevels = list(plot = rL_spatial_custom, sample = rL_nonspatial_custom)
   )

   nSamples <- 4
   fit_chain_2nf <- sampleMcmc(
      model, samples = nSamples, transient = 0, thin = 1,
      nChains = 1, adaptNf = c(0, 0), alignPost = FALSE, verbose = 0
   )
   fit_chain_3nf <- sampleMcmc(
      model, samples = nSamples, transient = 0, thin = 1,
      nChains = 1, adaptNf = c(0, 0), alignPost = FALSE, verbose = 0,
      initPar = list(Delta = list(matrix(1, 3, 1), matrix(1, 3, 1)))
   )

   for (ordering in list(
      list(fit = alignPosterior(c(fit_chain_2nf, fit_chain_3nf)), padded_chain = 1L),
      list(fit = alignPosterior(c(fit_chain_3nf, fit_chain_2nf)), padded_chain = 2L)
   )) {
      combined_fit <- ordering$fit
      padded_chain <- ordering$padded_chain
      for (posterior_sample in combined_fit$postList[[padded_chain]]) {
         expect_equal(length(posterior_sample$alphaInd[[1]]), 3L)
         expect_equal(length(posterior_sample$alpha[[1]]), 3L)
         expect_null(dim(posterior_sample$alphaInd[[1]]))
         expect_null(dim(posterior_sample$alpha[[1]]))
         expect_equal(posterior_sample$alphaInd[[1]][3], 1)
         expect_equal(posterior_sample$alpha[[1]][3], 0.1)

         expect_equal(length(posterior_sample$alphaInd[[2]]), 3L)
         expect_null(dim(posterior_sample$alphaInd[[2]]))
         expect_equal(posterior_sample$alphaInd[[2]][3], 1)
         expect_null(posterior_sample$alpha[[2]])
      }

      check_fitted_postprocessing(combined_fit, nSamples = nSamples, nChains = 2)
      coda_post <- convertToCodaObject(combined_fit)
      expect_true(all(as.matrix(coda_post$alpha[[1]][[padded_chain]])[, 3] == 0.1))
   }
})
