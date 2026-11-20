# The k + 1 largest of n observations from a generalized Pareto distribution
# with index gamma (survival (1 + gamma x)^(-1/gamma)), in decreasing order,
# by the Renyi representation of uniform order statistics.
top_gpd <- function(n, k, gamma) {
  e <- cumsum(stats::rexp(k + 1))
  p <- e / (e[k + 1] + stats::rgamma(1, n - k))
  if (gamma == 0) -log(p) else (p^(-gamma) - 1) / gamma
}
