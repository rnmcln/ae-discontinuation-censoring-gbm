## Scenario definitions
base_par <- list(
  pX = 0.4, bX = log(1.6), bU = 0.4, k = 1.3,
  gX = log(1.5), gU = 0.5, window = 12,
  aX = 0, aU = 0, fu_min = 12, fu_max = 36,
  nT = 60, nC = 180, median_C = 9
)

make_scenarios <- function() {
  main <- expand.grid(setting = "recurrent", median_C = 9, pD = c(0.10, 0.20, 0.30),
                      psi = c(1, 1.5, 2, 3), pi_bar = c(0.5, 1), aU = 0,
                      nT = 60, HR = c(1, 0.75), block = "main",
                      stringsAsFactors = FALSE)
  mnar <- expand.grid(setting = "recurrent", median_C = 9, pD = 0.20, psi = 2,
                      pi_bar = 0.5, aU = c(1, 2), nT = 60, HR = c(1, 0.75),
                      block = "mnar", stringsAsFactors = FALSE)
  size <- expand.grid(setting = "recurrent", median_C = 9, pD = 0.20, psi = 2,
                      pi_bar = 0.5, aU = 0, nT = c(40, 100), HR = c(1, 0.75),
                      block = "size", stringsAsFactors = FALSE)
  nd <- expand.grid(setting = "newly diagnosed", median_C = 16, pD = 0.20, psi = 2,
                    pi_bar = c(0.5, 1), aU = 0, nT = 60, HR = c(1, 0.75),
                    block = "newly_diagnosed", stringsAsFactors = FALSE)
  s <- rbind(main, mnar, size, nd)
  s$nC <- 3 * s$nT
  s$sid <- seq_len(nrow(s))
  s
}
