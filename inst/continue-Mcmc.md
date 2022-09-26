## Continuing MCMC in Hmsc

This document is temporary and discusses in detail the current
progress in enabling continuing MCMC chains in `sampleMcmcm`. Probably
this document will be removed once the functionality is finished.

The basic interface was already in `sampleMcmc` which had argument
`initPar` which is a named list of parameter values used to initialize
parameters instead of starting those parameters from priors. In the
old version, only one set of `initPar` could be supplied, and all
chains were started from the same `initPar`. The purpose of this
branch is to enable different `initPar` for each chain. Using the
values of last sample would -- in principle -- allow continuing
chains. However, there proved to be some "subtleties" making this a
bit tricky, like explained below in this document. Moreover, the goal
was to continue the *same* chain: that is, the continued chain should
be identical (within numerical precision) to having a longer sequence
originally. This brought along some other problems, but also allowed
testing the new functions: we should end up the same model when
running the same length sequence with the same RNG state in one or two
batches.

The major changes in the functions are:

- `sampleChain` accepts a list of `initPar` for each chain. Old
  behaviour is retained, and if only one set is provided, all chains
  will be intialized from that.
- New function `getLastPar` will extract the parameters of last sample
  in a form suitable for `initPar`. However, the last parameter values
  cannot be used directly as explained before, because various support
  functions modify the values returned with the `sampleMcmc`
  object. This is discussed in detail below in this document.
- New function `merge` appends continued chains to old
  chains. Function performs some sanity checks before merge, but
  usually it only warns of badly thought merges and makes them
  anyway. The only fatal error is that the number of chains do not
  match. Warnings are given if thin is changed, transient added or the
  `initPar` of new object are not equal to `getLastPar` parameters of
  old model. This may be changed, and some or all of these warnings
  may be changed to errors.

### Before sampling: `computeInitialParameters`

`sampleChain` will call `computeInitialParameters` before sampling
iterations. With new sampling, this will set initial values from
priors, but if `initPar` was supplied, the parameters in there will be
used as starting values. There were following problems:

- `Z` is not returned in `sampleMcmc` posterior, but it is calculated
  within function. `Z` is first calculated from fixed effects
  (`XScaled`, `Beta`) and random effects (`Lambda`, `Eta`) and then
  updated with `updateZ`. `updateZ` will call random numbers, and for
  determinisitic continuation of chain, RNG state must be set
  similarly as in the last iteration of old chain. **SOLVED:**
  `sampleChain` saves RNG state as attribute `Zstate` before last call
  to `updateZ` in sampling and in continued sampling sets RNG to this
  state before proper iterations and before calling
  `computeInitialParameters`.
  
- The previous approach will fail in RRR models with Poisson or
  log-Normal Poisson models. In these models `updateZ` needs also the
  `Z` of previous iteration, or the parameters of the iterations
  before the last one. **UNSOLVED:** We would need to access
  iterations before the last, and this is not even saved when
  thin>1. Other distributions can be continued, but `distr="poisson"`
  and `distr="lognormal poisson"` cannot be continued. Probably we
  need to warn with Poisson component models.

### Iterations: `sampleChain`

Predictable iterations must be started from the RNG state at the end
of saved sampling. **SOLVED:** The state of RNG is saved after the
last iteration as attribute `RNGstate` and RNG is set to this state
when continuing sampling.

### Storing samples: `combineParameters`

Sampling (in `sampleChain`) uses scaled data `XScaled`, but the
parameters are back-transformed to original non-scaled `X` in
`combineParameters` when storing the returned samples. Therefore
sampling cannot be continued from the values returned in the last
sample. 

- **SOLVED:** Function `getLastPar` will re-do the scaling of
  parameters `Beta`, `Gamma` and `V` back to internal forms, and these
  can be used to continue sampling.

- **NOTE:** RRR models are also handled. However, `combineParameters`
  does not back-transform RRR axes correctly: it uses `XRRRScalePar`
  which scales the original candidate variables of the RRR model, but
  the coefficients are for the latent RRR factors which have different
  mean and sd. Function `combineParameters` will use intercept
  coefficients (`m=0`, `s=1`) for the first RRR factor, and parameters
  of original variables for next axes. `getLastPar` ùses the same
  `XRRRScalePar` and correctly back-transform the coefficients back to
  internal form (even when this scaling is incorrect).

### After sampling: `alignPosterior`

Signs of latent random vectors `Lambda` and `Eta` are arbitrary. After
finishing sampling and exiting `sampleChain`, function
`alignPosterior` may be called (default is to call five times in a
loop) to align signs among samples. Switching signs similarly in
`Lambda` and `Eta` vectors shall not influence random effects, but it
_does_ influence updaters (I think it shouldn't buit it does).

- **SOLVED:** Information on sign switches is saved in attribute
  `alignment` for each chain, and this attribute is used to reset the
  signs similarly as in the original sampling. This provides
  predictable sequence. Currently the attribute saves the sign
  switches for each sample although only the last one is needed for
  continuing chains. This was done as I was curious to see how often
  vectors flip in sampling, but can be changed later.
  
- **NOT YET STUDIED:** Chains can have different numbers of latent
  vectors in `Lambda` and `Eta`, and these are padded to same
  dimensions. I have not yet studied how this influences getting and
  setting `initPar`.
  
- **NOT YET IMPLEMENTED:** in RRR models `alignPosterior` will flip
  latent weight vectors `wRRR` and parameters `Beta`, `Gamma`, and
  `V`. However, in my first tests these change rarely, and I haven't
  yet implemented this.
