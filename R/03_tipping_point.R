## Delta-adjusted multiple imputation (tipping-point analysis) for scenarios in
## which follow-up ends at every AE-related discontinuation (pi = 1).
## Usage: Rscript R/03_tipping_point.R [n_rep] [n_cores] [scenario ids, comma-separated]
source("R/functions.R"); source("R/scenarios.R")
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores() - 1L)
deltas <- c(1, 1.5, 2, 2.5, 3, 4, 5)
sc <- read.csv("results/scenarios_truth.csv")
sel <- sc[sc$block == "main" & sc$pi_bar == 1 & sc$pD == 0.2, ]
only <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else NULL
if (!is.null(only)) sel <- sel[sel$sid %in% only, ]
out <- list()
for (i in seq_len(nrow(sel))) {
  s <- sel[i, ]
  p <- modifyList(base_par, as.list(s[, c("psi", "b", "lamD", "aU", "a0", "nT", "nC")]))
  p$bA <- log(s$HR)
  set.seed(9000L + s$sid)
  seeds <- sample.int(.Machine$integer.max, n_rep)
  res <- mclapply(seq_len(n_rep), function(r) {
    set.seed(seeds[r])
    d <- generate_trial(p)
    do.call(rbind, lapply(deltas, function(dl) {
      est <- tryCatch(delta_mi(d, dl, M = 10),
                      error = function(e) c(logHR = NA, se = NA, S12.est = NA))
      data.frame(delta = dl, t(est), rep = r)
    }))
  }, mc.cores = n_cores)
  res <- do.call(rbind, res); res$sid <- s$sid
  out[[i]] <- res
  cat(sprintf("tipping point: scenario %d done\n", s$sid))
}
new <- do.call(rbind, out)
f <- "results/raw/tipping_point.rds"
if (!is.null(only) && file.exists(f)) {
  old <- readRDS(f)
  new <- rbind(old[!old$sid %in% only, ], new)
}
saveRDS(new[order(new$sid, new$rep, new$delta), ], f)
