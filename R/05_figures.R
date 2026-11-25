## Figures. Writes vector PDF/EPS and 600 dpi TIFF/PNG to the directory given
## as the first argument (default: figures/). Run with a UTF-8 locale.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(survival)})
source("R/functions.R"); source("R/scenarios.R")
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pf <- read.csv("results/performance.csv")
sc <- read.csv("results/scenarios_truth.csv")
fam <- "Liberation Sans"
ink <- "#2E3440"; mid <- "#4C566A"; light <- "#B5BBC6"; faint <- "#D8DEE9"; accent <- "#5E81AC"
col <- c(oracle = light, naive = accent, ipcw_base = mid, ipcw_disc = ink)
shp <- c(oracle = 1, naive = 16, ipcw_base = 17, ipcw_disc = 0)
lab <- c(oracle = "Complete follow-up", naive = "Naive",
         ipcw_base = "Baseline-weighted IPCW", ipcw_disc = "Discontinuation-weighted IPCW")
th <- theme_classic(base_size = 7, base_family = fam) +
  theme(axis.line = element_line(linewidth = 0.25, colour = "#3B4252"),
        axis.ticks = element_line(linewidth = 0.25, colour = "#3B4252"),
        axis.ticks.length = unit(1.2, "pt"),
        axis.text = element_text(colour = "#3B4252", size = 6),
        axis.title = element_text(colour = ink, size = 6.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, colour = ink, hjust = 0),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 6.5, colour = ink),
        panel.spacing.x = unit(10, "pt"), panel.spacing.y = unit(8, "pt"),
        legend.position = "none", plot.margin = margin(4, 70, 4, 4),
        plot.tag = element_text(size = 8, face = "plain", colour = ink))
save_fig <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(outdir, paste0(name, ".eps")), p, width = w_mm, height = h_mm, units = "mm",
         device = cairo_ps, fallback_resolution = 600)
  ggsave(file.path(outdir, paste0(name, ".tiff")), p, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, compression = "lzw")
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 600)
}
spread <- function(y, gap) {
  o <- order(y); ys <- y[o]
  for (i in seq_along(ys)[-1]) if (ys[i] - ys[i - 1] < gap) ys[i] <- ys[i - 1] + gap
  y[o] <- ys; y
}
xsc <- scale_x_continuous(breaks = c(1, 1.5, 2, 3), labels = c("1", "1.5", "2", "3"),
                          expand = expansion(add = c(0.15, 0.15)))
dodge <- position_dodge(width = 0.28)
main <- subset(pf, block == "main")
main$pD_lab <- factor(sprintf("Discontinuation %d%%", round(100 * main$pD)),
                      levels = sprintf("Discontinuation %d%%", c(10, 20, 30)))
main$pi_lab <- factor(ifelse(main$pi_bar == 1, "Follow-up ends at discontinuation", "Half of discontinuers lost"),
                      levels = c("Half of discontinuers lost", "Follow-up ends at discontinuation"))
main$method <- factor(main$method, levels = names(lab))

grid_panel <- function(d, yvar, lo, hi, href, ylab, gap, yfmt = waiver(), rows = "pi_lab") {
  d$y <- d[[yvar]]; d$lo <- lo; d$hi <- hi
  e <- d[d$psi == 3 & d$pD == 0.3 & is.finite(d$y), ]
  e <- do.call(rbind, lapply(split(e, e[[rows]], drop = TRUE), function(z) { z$ylab <- spread(z$y, gap); z }))
  ggplot(d, aes(psi, y, colour = method, shape = method, group = method)) +
    geom_hline(yintercept = href, linewidth = 0.25, colour = faint) +
    geom_line(linewidth = 0.35, position = dodge) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, position = dodge) +
    geom_point(size = 1.0, stroke = 0.4, position = dodge) +
    geom_text(data = e, aes(x = 3.28, y = ylab, label = lab[as.character(method)]),
              hjust = 0, size = 2.0, family = fam) +
    facet_grid(as.formula(paste(rows, "~ pD_lab")), switch = "y") + xsc +
    scale_y_continuous(labels = yfmt) +
    scale_colour_manual(values = col) + scale_shape_manual(values = shp) +
    coord_cartesian(clip = "off") +
    labs(x = "Hazard ratio for death after AE-related discontinuation (ψ)", y = ylab) + th
}

## Figure 1: mechanism, simulated patients and Kaplan-Meier curves (null scenario,
## 20% discontinuation, psi = 2, follow-up ends at discontinuation)
s20 <- sc[sc$sid == 20, ]
p20 <- modifyList(base_par, as.list(s20[, c("psi", "b", "lamD", "aU", "a0")]))
p20$bA <- 0
set.seed(2026)
tr <- observe(simulate_arm(400, 1L, p20), p20, lose = TRUE)
ec <- observe(simulate_arm(400, 0L, p20), p20, lose = FALSE)
pick <- function(d) {
  a <- d[d$disc == 1, ][1:6, ]; b <- d[d$disc == 0, ][1:6, ]
  x <- rbind(a, b); x <- x[order(pmin(x$T, x$C)), ]; x$row <- seq_len(nrow(x)); x
}
sw <- rbind(transform(pick(tr), cohort = "Single-arm trial"), transform(pick(ec), cohort = "External control"))
sw$cohort <- factor(sw$cohort, levels = c("Single-arm trial", "External control"))
sw$end_full <- pmin(sw$T, sw$C); sw$dead <- sw$T <= sw$C
unobs <- sw[sw$lost == 1, ]
pA <- ggplot(sw) +
  geom_segment(aes(x = 0, xend = time, y = row, yend = row), linewidth = 0.45, colour = mid) +
  geom_segment(data = unobs, aes(x = Dobs, xend = end_full, y = row, yend = row),
               linewidth = 0.35, colour = light, linetype = "22") +
  geom_point(data = subset(sw, disc == 1), aes(x = Dobs, y = row), shape = 21, fill = "white",
             colour = accent, size = 1.4, stroke = 0.5) +
  geom_point(data = subset(sw, dead & lost == 0), aes(x = T, y = row), shape = 4, colour = ink, size = 1.1, stroke = 0.45) +
  geom_point(data = subset(unobs, dead), aes(x = T, y = row), shape = 4, colour = light, size = 1.1, stroke = 0.45) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 6), expand = expansion(add = c(0, 1))) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Months from enrolment", y = "Simulated patients") + th +
  theme(axis.line.y = element_blank(), plot.margin = margin(4, 8, 4, 4))
## direct annotation in the trial panel
ann <- data.frame(cohort = factor("Single-arm trial", levels = levels(sw$cohort)),
                  x = c(12, 12, 12), y = c(3.2, 2.1, 1.0),
                  t = c("open circle: AE-related discontinuation", "x: death", "dashed grey: not observed after loss"))
pA <- pA + geom_text(data = ann, aes(x = x, y = y, label = t), hjust = 0, size = 1.9, family = fam, colour = mid)
set.seed(2027)
trL <- observe(simulate_arm(5000, 1L, p20), p20, lose = TRUE)
ecL <- observe(simulate_arm(5000, 0L, p20), p20, lose = FALSE)
km_df <- function(time, status, who) {
  f <- survfit(Surv(time, status) ~ 1)
  data.frame(t = c(0, f$time), s = c(1, f$surv), who = who)
}
kd <- rbind(km_df(trL$time, trL$status, "Trial, naive"),
            km_df(trL$time_full, trL$status_full, "Trial, complete follow-up"),
            km_df(ecL$time, ecL$status, "External control"))
kd <- kd[kd$t <= 30, ]
kcol <- c("Trial, naive" = accent, "Trial, complete follow-up" = light, "External control" = ink)
klty <- c("Trial, naive" = "solid", "Trial, complete follow-up" = "22", "External control" = "solid")
endp <- do.call(rbind, lapply(split(kd, kd$who), function(z) z[which.max(z$t), ]))
endp$y <- spread(endp$s, 0.05)
pC <- ggplot(kd, aes(t, s, colour = who, linetype = who)) +
  geom_step(linewidth = 0.45) +
  geom_text(data = endp, aes(x = 30.6, y = y, label = who), hjust = 0, size = 2.0, family = fam) +
  scale_colour_manual(values = kcol) + scale_linetype_manual(values = klty) +
  scale_x_continuous(breaks = seq(0, 30, 6), expand = expansion(add = c(0, 0.5))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) paste0(100 * x, "%")) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 1), clip = "off") +
  labs(x = "Months from enrolment", y = "Overall survival") + th + theme(plot.margin = margin(4, 95, 4, 4))
p1 <- (pA / pC) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "A")
save_fig(p1, "Fig1", 174, 135)

## Figure 2: hazard ratio and type I error under the null
d2 <- subset(main, HR == 1)
gm <- exp(d2$logHR_mean)
pH <- grid_panel(transform(d2, gm_tmp = gm), "gm_tmp", exp(d2$logHR_mean - 1.96 * d2$logHR_bias_mcse),
                 exp(d2$logHR_mean + 1.96 * d2$logHR_bias_mcse), 1, "Geometric mean of estimated HR", 0.022)
pT <- grid_panel(d2, "reject", pmax(0, d2$reject - 1.96 * d2$reject_mcse), d2$reject + 1.96 * d2$reject_mcse,
                 0.05, "Type I error (nominal 5%)", 0.019, yfmt = function(x) paste0(round(100 * x), "%"))
p2 <- (pH / pT) + plot_annotation(tag_levels = "A")
save_fig(p2, "Fig2", 174, 190)

## Figure 3: clinical scales, follow-up ends at discontinuation (pi = 1), null
d3 <- subset(main, HR == 1 & pi_bar == 1 & method != "ipcw_disc")
mk <- function(v, lo, hi, href, ylab, gap, yfmt = waiver()) {
  d <- d3; d$v <- v; d$lo <- lo; d$hi <- hi
  e <- d[d$psi == 3 & d$pD == 0.3, ]; e$ylab <- spread(e$v, gap)
  ggplot(d, aes(psi, v, colour = method, shape = method, group = method)) +
    geom_hline(yintercept = href, linewidth = 0.25, colour = faint) +
    geom_line(linewidth = 0.35, position = dodge) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, position = dodge) +
    geom_point(size = 1.0, stroke = 0.4, position = dodge) +
    geom_text(data = e, aes(x = 3.28, y = ylab, label = lab[as.character(method)]), hjust = 0, size = 2.0, family = fam) +
    facet_grid(. ~ pD_lab) + xsc + scale_y_continuous(labels = yfmt) +
    scale_colour_manual(values = col) + scale_shape_manual(values = shp) + coord_cartesian(clip = "off") +
    labs(x = "Hazard ratio for death after AE-related discontinuation (ψ)", y = ylab) + th
}
q1 <- mk(100 * d3$S12_bias, 100 * d3$S12_bias - 196 * d3$S12_bias_mcse, 100 * d3$S12_bias + 196 * d3$S12_bias_mcse,
         0, "Bias in 12-month OS (percentage points)", 0.9)
q2 <- mk(d3$median_bias, d3$median_bias - 1.96 * d3$median_bias_mcse, d3$median_bias + 1.96 * d3$median_bias_mcse,
         0, "Bias in median OS (months)", 0.17)
q3 <- mk(d3$S12_benchmark_success, pmax(0, d3$S12_benchmark_success - 1.96 * d3$S12_benchmark_mcse),
         d3$S12_benchmark_success + 1.96 * d3$S12_benchmark_mcse, 0.025,
         "Benchmark comparison declared positive", 0.012, yfmt = function(x) paste0(round(100 * x), "%"))
p3 <- (q1 / q2 / q3) + plot_annotation(tag_levels = "A")
save_fig(p3, "Fig3", 174, 175)

## Figure 4: tipping-point analysis
if (file.exists("results/tipping_point_performance.csv")) {
  tp <- read.csv("results/tipping_point_performance.csv")
  tp$HR_mean <- exp(tp$logHR_mean); tp$HR_true <- exp(tp$logHR_true)
  tp$panel <- factor(ifelse(tp$HR == 1, "True conditional HR 1.00", "True conditional HR 0.75"),
                     levels = c("True conditional HR 1.00", "True conditional HR 0.75"))
  tp$lo <- exp(tp$logHR_mean - 1.96 * tp$logHR_bias_mcse); tp$hi <- exp(tp$logHR_mean + 1.96 * tp$logHR_bias_mcse)
  tcol <- c(`1` = light, `1.5` = "#7B88A1", `2` = accent, `3` = ink)
  tp$psi_f <- factor(tp$psi)
  truth <- unique(tp[, c("panel", "psi_f", "HR_true")])
  ends <- tp[tp$delta == max(tp$delta), ]; ends$lab <- paste0("ψ = ", ends$psi)
  p4 <- ggplot(tp, aes(delta, HR_mean, colour = psi_f, group = psi_f)) +
    geom_hline(data = truth, aes(yintercept = HR_true, colour = psi_f), linewidth = 0.25, linetype = "22") +
    geom_line(linewidth = 0.35) + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35) +
    geom_point(size = 0.8) +
    geom_point(data = subset(tp, delta == psi), shape = 21, fill = "white", size = 1.7, stroke = 0.45) +
    geom_text(data = ends, aes(x = delta + 0.15, label = lab), hjust = 0, size = 2.1, family = fam) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(breaks = c(1, 1.5, 2, 2.5, 3, 4, 5), labels = c("1", "1.5", "2", "2.5", "3", "4", "5")) +
    scale_colour_manual(values = tcol) + coord_cartesian(clip = "off") +
    labs(x = "Assumed post-loss hazard multiplier δ in delta-adjusted imputation",
         y = "Geometric mean of estimated HR") + th
  save_fig(p4, "Fig4", 174, 70)
}
