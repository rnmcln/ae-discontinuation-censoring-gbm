# Supplementary Fig. S1: record flow of the literature review, drawn from review/review_summary.csv
suppressPackageStartupMessages(library(ggplot2))
s <- read.csv("review/review_summary.csv", stringsAsFactors = FALSE)
g <- function(k) s$n[s$item == k]
ink <- "#2E3440"; mid <- "#4C566A"; line <- "#4C566A"; fam <- "Liberation Sans"
ex <- s[grepl("^excluded: ", s$item), ]
ex$stage <- sub("^excluded: (.*) \\| .*$", "\\1", ex$item); ex$reason <- sub("^.* \\| ", "", ex$item)
lst <- function(st) { e <- ex[ex$stage == st, ]; e <- e[order(-e$n), ]; paste0(e$reason, " (", e$n, ")", collapse = "\n") }
tot <- s$N[s$item == "records: Full text | Included"][1]
n_type <- g("records: Article type | Excluded"); n_ta <- g("records: Title and abstract | Excluded"); n_ft <- g("records: Full text | Excluded")
n_inc <- g("records: Full text | Included"); n_comp <- g("records: Full text | Included as companion report")
main <- data.frame(y = c(10, 7.4, 4.8, 2.0),
  lab = c(paste0("Records identified in PubMed Central\n(search of 25 November 2026)\nn = ", tot),
          paste0("Records screened by title and abstract\nn = ", tot - n_type),
          paste0("Full texts assessed\nn = ", tot - n_type - n_ta),
          paste0("Included: ", n_inc + n_comp, " reports of ", n_inc, " trials\n(", n_comp, " companion follow-up report)")))
side <- data.frame(y = c(8.7, 6.1, 3.4),
  lab = c(paste0("Excluded by article type, n = ", n_type, "\n", lst("Article type")),
          paste0("Excluded, n = ", n_ta, "\n", lst("Title and abstract")),
          paste0("Excluded, n = ", n_ft, "\n", lst("Full text"))))
p <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 9.62, yend = 2.42), colour = line, linewidth = 0.3, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_segment(data = side, aes(x = 0, xend = 0.55, y = y, yend = y), colour = line, linewidth = 0.3, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_label(data = main, aes(x = 0, y = y, label = lab), family = fam, size = 2.3, colour = ink, label.size = 0.25, label.r = unit(0, "pt"), label.padding = unit(2.5, "pt"), lineheight = 1.05) +
  geom_label(data = side, aes(x = 0.58, y = y, label = lab), hjust = 0, family = fam, size = 2.1, colour = mid, label.size = 0.25, label.r = unit(0, "pt"), label.padding = unit(2.5, "pt"), lineheight = 1.05) +
  coord_cartesian(xlim = c(-0.85, 2.75), ylim = c(1.45, 10.55), clip = "off") +
  theme_void(base_family = fam) + theme(plot.margin = margin(4, 4, 4, 4))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/FigS1.pdf", p, width = 174, height = 125, units = "mm", device = cairo_pdf)
ggsave("figures/FigS1.png", p, width = 174, height = 125, units = "mm", dpi = 300, bg = "white")
ggsave("figures/FigS1.tiff", p, width = 174, height = 125, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("FigS1 written\n")
