# Quantile tables of the limiting null distributions of the tests of the
# extreme value condition, stored in R/sysdata.rda.
#
# * E (Dietrich, de Haan & Husler 2002, eq. 6; Husler & Li 2006), eta = 2,
#   for gamma_- on a grid in [-0.7, 0], with the sign of the J term as
#   corrected by Husler & Li (2006).
# * PE (Husler & Li 2006; T of Dietrich et al. 2002), eta = 0.5, 1, 2
#   (gamma >= 0; the limit does not depend on gamma).
# * kT_n (Drees, de Haan & Li 2006), eta = 1, maximum likelihood estimators,
#   for gamma on a grid in [-0.49, 2].
#
# Each limit is an integral of a function of a Brownian motion on (0, 1),
# approximated by the midpoint rule on m points. The generators reproduce
# Table 1 and 2 of Husler & Li (2006), Table 1 of Dietrich et al. (2002) and
# Table 1 of Drees, de Haan & Li (2006) (see exploracao/ and the reproduction
# vignette).
#
# Run from the package root (takes about half an hour):
#   Rscript data-raw/limit-tables.R

set.seed(20261120)
NSIM <- 50000L
M <- 4000L
PROBS <- c(seq(0.001, 0.999, by = 0.001), 0.9995)

sim_bm <- function(n, m) {
  inc <- matrix(stats::rnorm(m * n, sd = sqrt(1 / m)), m, n)
  w <- apply(inc, 2, cumsum)
  w0 <- rbind(0, w[-m, , drop = FALSE])
  list(mid = (w0 + w) / 2, w1 = w[m, ])
}

limit_pe <- function(n, m, eta) {
  tm <- (seq_len(m) - 0.5) / m
  b <- sim_bm(n, m)
  a <- b$mid / tm - rep(b$w1, each = m)
  ia <- colMeans(a)
  colMeans((a + outer(log(tm), ia))^2 * tm^eta)
}

limit_e <- function(n, m, g, eta) {
  tm <- (seq_len(m) - 0.5) / m
  b <- sim_bm(n, m)
  if (g == 0) {
    cfun <- -log(tm)
    jfun <- -(log(tm))^2 / 2
    aexp <- 1 / tm
  } else {
    cfun <- (tm^(-g) - 1) / g
    jfun <- tm^(-g) * log(tm) / g + (tm^(-g) - 1) / g^2
    aexp <- tm^(-g - 1)
  }
  a <- b$mid * aexp - rep(b$w1, each = m)
  p <- colMeans(a)
  q <- 2 * colMeans(cfun * a)
  r <- (1 - g)^2 * (1 - 2 * g) * ((0.5 - g) * q - 2 * p)
  colMeans(((1 - g) * a - (1 - g)^2 * outer(cfun, p) + outer(cfun, r) +
              (1 - g) * outer(jfun, r))^2 * tm^eta)
}

limit_t <- function(n, m, g, eta = 1) {
  tm <- (seq_len(m) - 0.5) / m
  b <- sim_bm(n, m)
  w <- b$mid
  w1 <- b$w1
  r <- colMeans(w / tm)
  if (abs(g) < 1e-12) {
    gam <- -colMeans((2 + log(tm)) * w / tm) + w1
    alp <- colMeans((3 + log(tm)) * w / tm) - 2 * w1
    l <- outer(tm, -w1) + outer(-0.5 * log(tm)^2 * tm, gam) +
      outer(log(tm) * tm, alp)
  } else {
    s <- colMeans(tm^(g - 1) * w)
    gam <- -((g + 1)^2 / g) * ((2 * g + 1) * s - r) + (g + 1) * w1
    alp <- -((g + 1) / g) * (r - (g + 1) * (2 * g + 1) * s) - (g + 2) * w1
    l <- outer(tm / g, gam / g - alp) + outer(tm * log(tm) / g, gam) -
      outer(tm^(1 + g) / g, g * w1 + gam / g - alp)
  }
  colMeans((w + l)^2 * tm^(eta - 2))
}

chunked <- function(fun, ..., m = M) {
  chunk <- max(50L, floor(2e6 / m))
  out <- numeric(0)
  while (length(out) < NSIM) {
    out <- c(out, fun(min(chunk, NSIM - length(out)), m, ...))
  }
  stats::quantile(out, PROBS, names = FALSE)
}

GAMMA_E <- seq(0, -0.7, by = -0.05)
GAMMA_T <- c(2, 1.5, 1, 0.75, 0.5, 0.25, 0.1, 0, -0.1, -0.2, -0.25, -0.3,
             -0.375, -0.45, -0.49)
ETA_PE <- c(0.5, 1, 2)

t0 <- Sys.time()
tab_e <- vapply(GAMMA_E, function(g) chunked(limit_e, g = g, eta = 2),
                numeric(length(PROBS)))
cat("E done", format(Sys.time() - t0), "\n")
# for eta <= 1 the integrand of PE (with log t) needs a finer grid near 0
tab_pe <- vapply(ETA_PE, function(eta) {
  chunked(limit_pe, eta = eta, m = if (eta <= 1) 20000L else M)
}, numeric(length(PROBS)))
cat("PE done", format(Sys.time() - t0), "\n")
# the limit of k T has log terms near 0 as well: finer grid
tab_t <- vapply(GAMMA_T, function(g) chunked(limit_t, g = g, m = 20000L),
                numeric(length(PROBS)))
cat("T done", format(Sys.time() - t0), "\n")

limit_tables <- list(
  probs = PROBS,
  E = list(gamma = GAMMA_E, q = tab_e, eta = 2),
  PE = list(eta = ETA_PE, q = tab_pe),
  T = list(gamma = GAMMA_T, q = tab_t, eta = 1),
  nsim = NSIM, m = M
)
usethis::use_data(limit_tables, internal = TRUE, overwrite = TRUE)
