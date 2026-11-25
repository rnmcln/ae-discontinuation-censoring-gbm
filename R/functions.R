## Core functions: data-generating mechanism, calibration, estimators, performance.
## All time units are months.

suppressPackageStartupMessages(library(survival))

## ---------------------------------------------------------------------------
## Data-generating mechanism
## ---------------------------------------------------------------------------
## Death hazard:  h(t) = h0(t) * exp(bA*A + bX*X + bU*U) * psi^{1[t > D]}
##   h0 Weibull with shape k and scale b: H0(t) = (t / b)^k
## AE-related treatment discontinuation (both arms, within a 12-month window):
##   hazard lamD * exp(gX*X + gU*U); D = Inf if it would occur after 12 months
## Loss of trial follow-up at discontinuation (single-arm trial only):
##   P(lost | D observed) = plogis(a0 + aX*X + aU*U)
## Administrative censoring: uniform accrual over 24 months, analysis 12 months
##   after last enrolment, i.e. C ~ U(12, 36), identical in both data sources.

death_time <- function(E, lp, D, psi, k, b) {
  H_D <- ifelse(is.finite(D), (D / b)^k, Inf)
  cH  <- exp(lp)
  pre <- E <= cH * H_D
  t_pre  <- b * (E / cH)^(1 / k)
  t_post <- b * (H_D + (E - cH * H_D) / (cH * psi))^(1 / k)
  ifelse(pre, t_pre, t_post)
}

simulate_arm <- function(n, A, p) {
  X <- rbinom(n, 1, p$pX)
  U <- rnorm(n)
  rateD <- p$lamD * exp(p$gX * X + p$gU * U)
  D <- rexp(n, rateD)
  D[D > p$window] <- Inf
  E <- rexp(n)
  lp <- p$bA * A + p$bX * X + p$bU * U
  Tt <- death_time(E, lp, D, p$psi, p$k, p$b)
  C  <- runif(n, p$fu_min, p$fu_max)
  data.frame(A = A, X = X, U = U, D = D, T = Tt, C = C)
}

## Apply observation scheme. `lose` = TRUE for the single-arm trial.
observe <- function(d, p, lose) {
  disc <- is.finite(d$D) & d$D < pmin(d$T, d$C)
  lost <- rep(FALSE, nrow(d))
  if (lose) {
    pl <- plogis(p$a0 + p$aX * d$X + p$aU * d$U)
    lost <- disc & (runif(nrow(d)) < pl)
  }
  t_full <- pmin(d$T, d$C)
  s_full <- as.integer(d$T <= d$C)
  d$disc   <- as.integer(disc)
  d$lost   <- as.integer(lost)
  d$Dobs   <- ifelse(disc, d$D, NA_real_)
  d$time_full   <- t_full            # complete ascertainment (oracle)
  d$status_full <- s_full
  d$time   <- ifelse(lost, d$D, t_full)
  d$status <- ifelse(lost, 0L, s_full)
  d
}

generate_trial <- function(p) {
  tr <- observe(simulate_arm(p$nT, 1L, p), p, lose = TRUE)
  ec <- observe(simulate_arm(p$nC, 0L, p), p, lose = FALSE)
  d <- rbind(tr, ec)
  d$id <- seq_len(nrow(d))
  d
}

## ---------------------------------------------------------------------------
## Calibration (common random numbers, large sample)
## ---------------------------------------------------------------------------
## Solve for Weibull scale b (control median OS = target) and discontinuation
## rate lamD (cumulative incidence of AE-related discontinuation by 12 months in
## the control population, death as competing event = target).
calibrate <- function(p, target_median, target_cif, n = 4e5, seed = 1) {
  set.seed(seed)
  X <- rbinom(n, 1, p$pX); U <- rnorm(n); Z <- rexp(n); E <- rexp(n)
  lp <- p$bX * X + p$bU * U
  sD <- exp(p$gX * X + p$gU * U)
  b <- target_median; lamD <- -log(1 - target_cif) / 12
  for (it in 1:30) {
    D <- Z / (lamD * sD); D[D > p$window] <- Inf
    f_med <- function(lb) median(death_time(E, lp, D, p$psi, p$k, exp(lb))) - target_median
    b <- exp(uniroot(f_med, log(c(0.2, 20) * target_median), tol = 1e-7)$root)
    f_cif <- function(ll) {
      Dl <- Z / (exp(ll) * sD); Dl[Dl > p$window] <- Inf
      Tl <- death_time(E, lp, Dl, p$psi, p$k, b)
      mean(Dl < pmin(Tl, 12)) - target_cif
    }
    lamD_new <- exp(uniroot(f_cif, log(c(1e-5, 2)), tol = 1e-8)$root)
    if (abs(lamD_new - lamD) / lamD < 1e-6) { lamD <- lamD_new; break }
    lamD <- lamD_new
  }
  list(b = b, lamD = lamD)
}

## Intercept a0 so that mean loss probability among discontinuers equals pi_bar
calibrate_a0 <- function(p, pi_bar, n = 4e5, seed = 2) {
  if (pi_bar >= 1) return(Inf)
  if (p$aX == 0 && p$aU == 0) return(qlogis(pi_bar))
  set.seed(seed)
  d <- observe(simulate_arm(n, 1L, modifyList(p, list(a0 = 0))), p, lose = FALSE)
  dd <- d[d$disc == 1, ]
  f <- function(a0) mean(plogis(a0 + p$aX * dd$X + p$aU * dd$U)) - pi_bar
  uniroot(f, c(-15, 15), tol = 1e-8)$root
}

## True treatment-policy estimands under complete ascertainment.
## Common random numbers: the same covariates, frailties and uniform draws are
## used for A = 1 and A = 0, so Monte Carlo error largely cancels in contrasts.
true_values <- function(p, n = 5e5, seed = 3) {
  set.seed(seed)
  X <- rbinom(n, 1, p$pX); U <- rnorm(n); Z <- rexp(n); E <- rexp(n)
  C <- runif(n, p$fu_min, p$fu_max)
  D <- Z / (p$lamD * exp(p$gX * X + p$gU * U)); D[D > p$window] <- Inf
  T1 <- death_time(E, p$bA + p$bX * X + p$bU * U, D, p$psi, p$k, p$b)
  T0 <- death_time(E, p$bX * X + p$bU * U, D, p$psi, p$k, p$b)
  d <- data.frame(time = c(pmin(T1, C), pmin(T0, C)),
                  status = c(as.integer(T1 <= C), as.integer(T0 <= C)),
                  A = rep(1:0, each = n))
  lhr <- if (p$bA == 0) 0 else unname(coef(coxph(Surv(time, status) ~ A, data = d, ties = "breslow")))
  c(logHR = lhr, S12_T = mean(T1 > 12), S12_C = mean(T0 > 12),
    cif_disc12_T = mean(is.finite(D) & D < pmin(T1, 12)),
    cif_disc12_C = mean(is.finite(D) & D < pmin(T0, 12)),
    med_T = median(T1), med_C = median(T0))
}

## ---------------------------------------------------------------------------
## Estimators
## ---------------------------------------------------------------------------
km_at <- function(fit, t) {
  s <- summary(fit, times = t, extend = TRUE)
  c(est = s$surv, se = s$std.err)
}

km_median <- function(fit) {
  q <- tryCatch(quantile(fit, probs = 0.5, conf.int = FALSE), error = function(e) NA_real_)
  as.numeric(q)[1]
}

cox_out <- function(fit) {
  b <- unname(coef(fit)); se <- sqrt(diag(vcov(fit)))[1]
  c(logHR = b, se = unname(se))
}

## Counting-process expansion of treated arm at a set of cut times
split_at <- function(d, cuts) {
  if (!length(cuts)) { d$tstart <- 0; d$tstop <- d$time; d$event <- d$status; return(d) }
  s <- survSplit(Surv(time, status) ~ ., data = d, cut = sort(unique(cuts)),
                 start = "tstart", end = "tstop", event = "event")
  s
}

analyse <- function(d, tau = 12) {
  out <- list()
  ## M0 oracle: complete follow-up of the trial arm
  f0 <- coxph(Surv(time_full, status_full) ~ A, data = d)
  k0 <- survfit(Surv(time_full, status_full) ~ 1, data = d[d$A == 1, ])
  out$oracle <- c(cox_out(f0), S12 = km_at(k0, tau), med = km_median(k0))
  ## M1 naive: loss at discontinuation treated as non-informative censoring
  f1 <- coxph(Surv(time, status) ~ A, data = d)
  k1 <- survfit(Surv(time, status) ~ 1, data = d[d$A == 1, ])
  out$naive <- c(cox_out(f1), S12 = km_at(k1, tau), med = km_median(k1))

  tr <- d[d$A == 1, ]; ec <- d[d$A == 0, ]
  n_lost <- sum(tr$lost); n_disc <- sum(tr$disc)
  ## M2 IPCW with baseline covariate: KM of the loss process stratified by X
  if (n_lost > 0) {
    cuts <- sort(unique(tr$time[tr$lost == 1]))
    st <- split_at(tr, cuts)
    st$w <- 1
    for (x in 0:1) {
      ix <- tr$X == x
      if (!any(ix)) next
      g <- survfit(Surv(time, lost) ~ 1, data = tr[ix, ])
      G <- stepfun(g$time, c(1, g$surv))
      jx <- st$X == x
      st$w[jx] <- 1 / pmax(G(st$tstart[jx]), 1e-3)
    }
    ec2 <- ec; ec2$tstart <- 0; ec2$tstop <- ec2$time; ec2$event <- ec2$status; ec2$w <- 1
    cols <- c("id", "A", "X", "tstart", "tstop", "event", "w")
    dd <- rbind(st[, cols], ec2[, cols])
    f2 <- coxph(Surv(tstart, tstop, event) ~ A, data = dd, weights = w, cluster = id)
    k2 <- survfit(Surv(tstart, tstop, event) ~ 1, data = dd[dd$A == 1, ], weights = w,
                  id = id, robust = TRUE)
    out$ipcw_base <- c(cox_out(f2), S12 = km_at(k2, tau), med = km_median(k2))
  } else out$ipcw_base <- out$naive

  ## M3 IPCW conditional on discontinuation: P(lost | discontinued, X)
  retained <- tr$disc == 1 & tr$lost == 0
  if (n_lost == 0) {
    out$ipcw_disc <- out$naive
  } else if (sum(retained) == 0) {
    out$ipcw_disc <- c(logHR = NA, se = NA, S12.est = NA, S12.se = NA, med = NA)
  } else {
    dd_disc <- tr[tr$disc == 1, ]
    m <- suppressWarnings(glm(lost ~ X, family = binomial, data = dd_disc))
    pl <- predict(m, newdata = tr, type = "response")
    cut_ids <- which(retained)
    st <- tr
    st$tstart <- 0; st$tstop <- st$time; st$event <- st$status; st$w <- 1
    a <- st[cut_ids, ]; bpart <- a
    a$tstop <- a$Dobs; a$event <- 0L
    bpart$tstart <- bpart$Dobs; bpart$w <- 1 / pmax(1 - pl[cut_ids], 1e-3)
    st <- rbind(st[-cut_ids, ], a, bpart)
    st <- st[st$tstop > st$tstart, ]
    ec2 <- ec; ec2$tstart <- 0; ec2$tstop <- ec2$time; ec2$event <- ec2$status; ec2$w <- 1
    cols <- c("id", "A", "X", "tstart", "tstop", "event", "w")
    dd <- rbind(st[, cols], ec2[, cols])
    f3 <- coxph(Surv(tstart, tstop, event) ~ A, data = dd, weights = w, cluster = id)
    k3 <- survfit(Surv(tstart, tstop, event) ~ 1, data = dd[dd$A == 1, ], weights = w,
                  id = id, robust = TRUE)
    out$ipcw_disc <- c(cox_out(f3), S12 = km_at(k3, tau), med = km_median(k3))
  }
  res <- do.call(rbind, out)
  data.frame(method = rownames(res), res, n_lost = n_lost, n_disc = n_disc,
             row.names = NULL, check.names = FALSE)
}

## Delta-adjusted multiple imputation for patients lost at discontinuation.
## Post-loss hazard = delta * hazard predicted by a Weibull PH model fitted to the
## observed trial-arm data (covariate X). M imputations, Rubin's rules on log HR.
delta_mi <- function(d, delta, M = 10, tau = 12) {
  tr <- d[d$A == 1, ]; ec <- d[d$A == 0, ]
  fit <- survreg(Surv(time, status) ~ X, data = tr, dist = "weibull")
  theta <- c(coef(fit), log(fit$scale))
  V <- vcov(fit)                                       # includes Log(scale)
  R <- chol(V)
  Xm <- cbind(1, tr$X)
  L <- which(tr$lost == 1)
  est <- se2 <- s12 <- numeric(M)
  for (m in seq_len(M)) {
    ## Proper imputation: draw parameters from their approximate sampling distribution
    th <- theta + drop(rnorm(length(theta)) %*% R)
    k <- 1 / exp(th[3])
    bsc <- exp(drop(Xm %*% th[1:2]))                   # S(t) = exp(-(t / bsc)^k)
    ti <- tr
    if (length(L)) {
      HD <- (tr$Dobs[L] / bsc[L])^k
      Tn <- bsc[L] * (HD + rexp(length(L)) / delta)^(1 / k)
      ti$time[L]   <- pmin(Tn, tr$C[L])
      ti$status[L] <- as.integer(Tn <= tr$C[L])
    }
    dm <- rbind(ti, ec)
    f <- coxph(Surv(time, status) ~ A, data = dm)
    est[m] <- coef(f); se2[m] <- vcov(f)[1, 1]
    s12[m] <- km_at(survfit(Surv(time, status) ~ 1, data = ti), tau)[1]
  }
  qbar <- mean(est); W <- mean(se2); B <- var(est)
  c(logHR = qbar, se = sqrt(W + (1 + 1 / M) * B), S12.est = mean(s12))
}

## ---------------------------------------------------------------------------
## Performance measures (Morris, White & Crowther 2019) with Monte Carlo SEs
## ---------------------------------------------------------------------------
perf <- function(est, se, truth, null_value = 0) {
  ok <- is.finite(est) & is.finite(se)
  e <- est[ok]; s <- se[ok]; n <- length(e)
  if (n < 2) return(c(n = n, bias = NA, bias_mcse = NA, empse = NA, modse = NA,
                      cover = NA, cover_mcse = NA, reject = NA, reject_mcse = NA))
  bias <- mean(e) - truth; empse <- sd(e)
  lo <- e - qnorm(0.975) * s; hi <- e + qnorm(0.975) * s
  cov <- mean(lo <= truth & truth <= hi)
  rej <- mean(abs((e - null_value) / s) > qnorm(0.975))
  c(n = n, bias = bias, bias_mcse = empse / sqrt(n), empse = empse,
    modse = sqrt(mean(s^2)), cover = cov, cover_mcse = sqrt(cov * (1 - cov) / n),
    reject = rej, reject_mcse = sqrt(rej * (1 - rej) / n))
}
