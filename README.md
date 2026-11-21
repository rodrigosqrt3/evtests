# evtests

[![R-CMD-check](https://github.com/rodrigosqrt3/evtests/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/evtests/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/evtests/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/evtests)

Tests and diagnostics for the upper tail of a distribution, to be used before
an extreme value analysis: tests of the extreme value condition, tests of the
domain of attraction and of a finite end-point, estimators of the extreme
value index and of the end-point, goodness of fit of the generalized Pareto
distribution, and the choice of the number of upper order statistics.

Every method is implemented from its source paper and tested against the
paper's numerical results: worked examples are recomputed from the printed
data, and simulation tables are re-run with the published design. Where a
printed value differs from the value given by the paper's own formulas or
design, the difference is documented in `vignette("reproduction")`.

Throughout the package, `k` is the number of excesses over the `(k + 1)`-th
largest observation.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rodrigosqrt3/evtests")
```

## Methods

| Topic | Functions |
|---|---|
| Extreme value index | `hill1975()`, `pickands1975()`, `dekkers1989()`, `drees2004()` (generalized Pareto maximum likelihood, computed with `grimshaw1993()`) |
| Extreme value condition | `dietrich2002()` (statistics `E` and `PE`), `drees2006()` |
| Gumbel domain (`gamma = 0`) | `hasofer_wang1992()`, `segers_teugels2000()`, `neves2006()`, `neves_fraga_alves2007()`, `fraga_alves2016()`, `castillo2014()`, `castillo2014_cv()` |
| Finite end-point | `neves_pereira2010()`, `hall_wang1999()`, `smith1985()`, `fraga_alves_neves2014()` (with the upper bound of Fraga Alves et al. 2016) |
| Generalized Pareto fit | `grimshaw1993()`, `choulakian_stephens2001()` |
| Choice of `k` or of the threshold | `wang1995()`, `drees_kaufmann1998()`, `danielsson2001()`, `guillou_hall2001()`, `bader2018()` |

## Example

```r
library(evtests)
# Hill (1975), Table 1: the 30 largest cities of the United States in 1940
hill1975(us_cities_1940$population, k = c(5, 10, 29))

# Hall & Wang (1999): lower end-point of the 100 m times
hall_wang1999(olympic_100m$time, k = 11, alpha = 0.1, side = "lower")
```
