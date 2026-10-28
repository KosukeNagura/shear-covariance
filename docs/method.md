# Current method and scope

## Observable and signal model

The current demonstration has ten angular bins. Its data vector ordering is all ξ+ bins followed by all ξ− bins. The signal-only model is an EE shear spectrum with BB set to zero, transformed on the sphere with the spin-2 correlation functions. Catalog signal maps use the established map convention and TreeCorr settings recorded in the simulation output.

## Shape-noise model

The catalog simulation adds independent, isotropic Gaussian noise to the two shear components of each selected galaxy. The component variance is calibrated from weighted catalog ellipticity moments and checked against the existing convergence noise convention. It records signal-only, noise-only, and signal-plus-noise measurements separately; total correlations are measured from the summed galaxy shear values.

For grid pixel (p), the variance of the weighted mean of independent galaxy noise is

\[
v_p=\sigma^2_{\rm component}\,
\frac{\sum_{g\in p}w_g^2}{(\sum_{g\in p}w_g)^2}.
\]

The grid covariance code adds (v_p I_2) to the same-pixel signal covariance block and evaluates the resulting disconnected Gaussian/Wick contractions. The current code uses unordered grid pairs (i<j), fixed full-pair denominators, and the established quartic-counts kernel.

## Pair thinning

The down-sampled calculation retains each pair with probability (p=0.001), repeats the thinning 100 times, and corrects sampled pair-pair sums by their inclusion probabilities. It also applies the self-pair correction required when both factors refer to the same retained pair. The output mean estimates the full-grid pair sum; the reported SEM measures variation across thinning realizations, not the finite-simulation error.

## What the current comparison establishes

The catalog TreeCorr simulation covariance and the grid estimator covariance are compared as a diagnostic. Differences can include estimator, pixelization, sampling, and finite-realization effects. Do not interpret their difference as a pure test of the Wick formula without additional matched-estimator checks.

The current signal model is E-only. The code has not yet been generalized to arbitrary EE/BB/EB inputs or to the HSC Final-Year B-mode data-vector definition.
