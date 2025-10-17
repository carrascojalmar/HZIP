# HZIP
## Likelihood-based inference for joint modeling of correlated count and binary outcomes with  extra variability and zeros

### Authors
Lizandra C. Fabio,
Jalmar M. F. Carrasco,
Víctor H. Lachos,
Ming-Hui Chen.

### Abstract
This paper introduces an inference approach for jointly modeling correlated count and binary outcomes, effectively capturing the dependence structure between subjects, overdispersion, and excess zeros. The methodology is based on a hierarchical zero-inflated Poisson (HZIP) model. The likelihood function is derived analytically using Newton’s binomial expansion and is expressed as a sum of the product of marginal distributions, which requires only a single integral over a set of binary variables. This formulation enables the simultaneous modeling of zero inflation from the Bernoulli component while providing a more accurate evaluation of the HZIP model's parsimony. Furthermore, maximum likelihood (ML) estimation under this inference framework is computationally efficient, as the necessary integrals can be solved analytically. Monte Carlo simulations are conducted to evaluate the performance of the maximum likelihood estimators under various scenarios. The proposed HZIP model is applied to the salamander mating dataset, demonstrating superior performance compared to the \textcolor{red}{hurdle Poisson model}. Finally, a residual analysis using randomized quantile residuals is proposed to assess potential departures from the HZIP model.

### To install the package is directly from the github.com repository, by running

R (>= 3.5)

require(devtools)

devtools::install_github("carrascojalmar/HZIP")
