# HZIP 0.1.5

   * Fixed the marginal likelihood of the zero-inflation component (`dBerLGG`):
     the alternating sum over s was taken over all observed zeros, including
     those with k_j = 1, which produced biased estimates whenever clusters
     contained zeros. The sum now runs only over positions with y_j = 0 and
     k_j = 0, matching the expansion of pi^(z - k) in Theorem 1.
   * Slow model fits in the `rhzip()` examples are now wrapped in `\donttest{}`.
   
# HZIP 0.1.4

* Fixed a bug in `rhzip()` where the Bernoulli draw used to generate
  structural zeros had its probability inverted (`rbinom(m[i], 1, p)`
  instead of `rbinom(m[i], 1, 1 - p)`), which produced
  `P(Y=0) = (1-pi) + pi*exp(-mu)` instead of the correct
  `P(Y=0) = pi + (1-pi)*exp(-mu)`.
* Fixed `simulate.hzip()`, which implemented a hurdle-Poisson generative
  mechanism (rejection/retry sampling to avoid zero counts) instead of the
  correct ZIP/HZIP mixture mechanism (an ordinary, non-truncated Poisson
  draw that can itself be zero). This affected the parametric bootstrap
  used by `testZI()` and `testDisp()`.
* Fixed `rhzip()` to no longer silently drop the first column of `x1`
  and/or `x2`. Previously the first column of each matrix was always
  removed, assuming it was an intercept; this discarded a real covariate
  whenever `x1`/`x2` did not include an intercept as their first column.
  All supplied columns are now kept and renamed (`x1` -> `x1, ..., xp1`;
  `x2` -> `w1, ..., wp2`).
* Corrected the `hzip()` documentation, which described the two-part
  formula as `y ~ w_zero + ... | x_count + ...`; the actual (and now
  documented) convention is `y ~ x_zero + ... | w_count + ...`, matching
  `rhzip()`'s column-naming convention.
* `AICz`, `BICz`, `AIC_HZIP` and `BIC_HZIP` are temporarily disabled
  (commented out) in `hzip()`, `summary.hzip()` and `print.summary.hzip()`
  and are no longer computed, returned, or printed.
* In `print.summary.hzip()`, `Convergence` is now shown as an additional
  column in the main summary table (alongside `logLik`, `AIC`, `BIC`,
  `Iter`) instead of on a separate line.

# HZIP 0.1.3

* The class of fitted objects was renamed from `HZIP` to `hzip`; all S3
  methods (`print`, `summary`, `residuals`, `simulate`, `testDisp`, `testZI`)
  were updated accordingly.
* `rHZIP()` was renamed to `rhzip()`.
* `envelope.HZIP()` was renamed to `envelope()`.
* `testDisp()` and `testZI()` are now S3 generics with explicit arguments
  (`nsim`, `alternative`, `alpha`, `seed`, `plot`).
* The `salamanders` dataset now includes the variable `sample`; examples
  were updated to use it.
* Improved package startup message and documentation.

# HZIP 0.1.2

* Previous CRAN release.
