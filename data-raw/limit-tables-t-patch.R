# Recomputes only the table of k T (Drees, de Haan & Li 2006) of
# R/sysdata.rda with the grid used in data-raw/limit-tables.R (m = 20000),
# keeping the other tables. Run from the package root.
src <- readLines("data-raw/limit-tables.R")
stop_at <- grep("^t0 <- Sys.time", src) - 1L
eval(parse(text = src[seq_len(stop_at)]))
set.seed(20261122)
load("R/sysdata.rda")
limit_tables$T$q <- vapply(GAMMA_T, function(g) {
  chunked(limit_t, g = g, m = 20000L)
}, numeric(length(PROBS)))
usethis::use_data(limit_tables, internal = TRUE, overwrite = TRUE)
