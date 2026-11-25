## Calibrate nuisance parameters and compute true estimands for every scenario.
source("R/functions.R"); source("R/scenarios.R")
sc <- make_scenarios()
key <- unique(sc[, c("median_C", "pD", "psi")])
cal <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  p <- modifyList(base_par, list(psi = key$psi[i], bA = 0))
  cb <- calibrate(p, key$median_C[i], key$pD[i], seed = 101)
  data.frame(key[i, ], b = cb$b, lamD = cb$lamD)
}))
sc <- merge(sc, cal, by = c("median_C", "pD", "psi"))
sc <- sc[order(sc$sid), ]
sc$a0 <- NA_real_
truth <- vector("list", nrow(sc))
for (i in seq_len(nrow(sc))) {
  p <- modifyList(base_par, as.list(sc[i, c("psi", "b", "lamD", "aU", "nT", "nC")]))
  p$bA <- log(sc$HR[i])
  sc$a0[i] <- calibrate_a0(p, sc$pi_bar[i], seed = 202)
  p$a0 <- sc$a0[i]
  truth[[i]] <- true_values(p, n = 1e6, seed = 303)
  cat(sprintf("scenario %d done\n", sc$sid[i]))
}
sc <- cbind(sc, do.call(rbind, truth))
dir.create("results", showWarnings = FALSE)
write.csv(sc, "results/scenarios_truth.csv", row.names = FALSE)
