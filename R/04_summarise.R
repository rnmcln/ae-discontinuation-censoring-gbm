## Performance measures for all scenarios and methods.
source("R/functions.R")
sc <- read.csv("results/scenarios_truth.csv")
files <- sprintf("results/raw/sim_%03d.rds", sc$sid)
rows <- list()
for (i in seq_len(nrow(sc))) {
  if (!file.exists(files[i])) next
  r <- readRDS(files[i]); s <- sc[i, ]
  ## The discontinuation-conditional IPCW estimator is not identified when no
  ## discontinuing patient can remain under follow-up (pi = 1).
  if (!"med" %in% names(r)) r$med <- NA_real_
  if (s$pi_bar >= 1) r[r$method == "ipcw_disc", c("logHR", "se", "S12.est", "S12.se", "med")] <- NA
  for (m in unique(r$method)) {
    x <- r[r$method == m, ]
    ph <- perf(x$logHR, x$se, s$logHR, null_value = 0)
    ps <- perf(x$S12.est, x$S12.se, s$S12_T, null_value = s$S12_C)
    ## Single-arm benchmark comparison: lower 95% limit of 12-month OS above the
    ## true 12-month OS of the control population (one-sided 2.5% level)
    okb <- is.finite(x$S12.est) & is.finite(x$S12.se)
    bench <- mean(x$S12.est[okb] - qnorm(0.975) * x$S12.se[okb] > s$S12_C)
    okm <- is.finite(x$med)
    med_bias <- mean(x$med[okm]) - s$med_T
    rows[[length(rows) + 1]] <- data.frame(
      sid = s$sid, method = m, n_rep = nrow(x),
      n_estimable = unname(ph["n"]), prop_not_estimable = 1 - unname(ph["n"]) / nrow(x),
      logHR_true = s$logHR, HR_true = exp(s$logHR),
      logHR_mean = mean(x$logHR, na.rm = TRUE),
      logHR_bias = unname(ph["bias"]), logHR_bias_mcse = unname(ph["bias_mcse"]),
      logHR_empse = unname(ph["empse"]), logHR_modse = unname(ph["modse"]),
      logHR_cover = unname(ph["cover"]), logHR_cover_mcse = unname(ph["cover_mcse"]),
      reject = unname(ph["reject"]), reject_mcse = unname(ph["reject_mcse"]),
      S12_true = s$S12_T,
      S12_bias = unname(ps["bias"]), S12_bias_mcse = unname(ps["bias_mcse"]),
      S12_empse = unname(ps["empse"]),
      S12_cover = unname(ps["cover"]), S12_cover_mcse = unname(ps["cover_mcse"]),
      S12_benchmark_success = bench, S12_benchmark_mcse = sqrt(bench * (1 - bench) / sum(okb)),
      median_true = s$med_T, median_bias = med_bias,
      median_bias_mcse = sd(x$med[okm]) / sqrt(sum(okm)), median_n = sum(okm),
      mean_n_disc = mean(x$n_disc), mean_n_lost = mean(x$n_lost),
      prop_no_loss = mean(x$n_lost == 0))
  }
}
pf <- do.call(rbind, rows)
pf <- merge(sc[, c("sid", "block", "setting", "median_C", "pD", "psi", "pi_bar", "aU", "nT", "nC", "HR")],
            pf, by = "sid")
pf <- pf[order(pf$sid, match(pf$method, c("oracle", "naive", "ipcw_base", "ipcw_disc"))), ]
num <- vapply(pf, is.numeric, logical(1)); pf[num] <- lapply(pf[num], signif, digits = 8)
write.csv(pf, "results/performance.csv", row.names = FALSE)

if (file.exists("results/raw/tipping_point.rds")) {
  tp <- readRDS("results/raw/tipping_point.rds")
  truth <- sc[, c("sid", "psi", "HR", "logHR", "S12_T")]
  names(truth)[names(truth) == "logHR"] <- "logHR_true"
  tp <- merge(tp, truth, by = "sid")
  agg <- do.call(rbind, lapply(split(tp, list(tp$sid, tp$delta), drop = TRUE), function(x) {
    ph <- perf(x$logHR, x$se, x$logHR_true[1], 0)
    data.frame(sid = x$sid[1], psi = x$psi[1], HR = x$HR[1], delta = x$delta[1],
               n = sum(is.finite(x$logHR)), logHR_true = x$logHR_true[1],
               logHR_mean = mean(x$logHR, na.rm = TRUE), logHR_bias = unname(ph["bias"]),
               logHR_bias_mcse = unname(ph["bias_mcse"]), logHR_cover = unname(ph["cover"]),
               reject = unname(ph["reject"]), reject_mcse = unname(ph["reject_mcse"]),
               S12_bias = mean(x$S12.est, na.rm = TRUE) - x$S12_T[1])
  }))
  agg <- agg[order(agg$sid, agg$delta), ]
  num <- vapply(agg, is.numeric, logical(1)); agg[num] <- lapply(agg[num], signif, digits = 8)
  write.csv(agg, "results/tipping_point_performance.csv", row.names = FALSE)
}
