## Monte Carlo simulation: n_rep replicates per scenario.
## Usage: Rscript R/02_run_simulation.R [n_rep] [n_cores]
source("R/functions.R"); source("R/scenarios.R")
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 2000L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores() - 1L)
sc <- read.csv("results/scenarios_truth.csv")
dir.create("results/raw", recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(sc))) {
  f <- sprintf("results/raw/sim_%03d.rds", sc$sid[i])
  if (file.exists(f)) next
  p <- modifyList(base_par, as.list(sc[i, c("psi", "b", "lamD", "aU", "a0", "nT", "nC")]))
  p$bA <- log(sc$HR[i])
  set.seed(20251125L + sc$sid[i])
  seeds <- sample.int(.Machine$integer.max, n_rep)
  res <- mclapply(seq_len(n_rep), function(r) {
    set.seed(seeds[r])
    d <- generate_trial(p)
    out <- tryCatch(analyse(d), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    out$rep <- r
    out
  }, mc.cores = n_cores)
  res <- do.call(rbind, res)
  res$sid <- sc$sid[i]
  saveRDS(res, f)
  cat(sprintf("%s scenario %d/%d complete (%d rows)\n", format(Sys.time(), "%H:%M:%S"),
              i, nrow(sc), nrow(res)))
}
