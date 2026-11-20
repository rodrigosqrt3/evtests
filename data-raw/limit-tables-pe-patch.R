# Recomputes only the PE columns of R/sysdata.rda with the grid used in
# data-raw/limit-tables.R (m = 20000 for eta <= 1), keeping the other
# tables. Run from the package root.
src <- readLines("data-raw/limit-tables.R")
stop_at <- grep("^t0 <- Sys.time", src) - 1L
eval(parse(text = src[seq_len(stop_at)]))
set.seed(20261121)
load("R/sysdata.rda")
limit_tables$PE$q <- vapply(ETA_PE, function(eta) {
  chunked(limit_pe, eta = eta, m = if (eta <= 1) 20000L else M)
}, numeric(length(PROBS)))
usethis::use_data(limit_tables, internal = TRUE, overwrite = TRUE)
