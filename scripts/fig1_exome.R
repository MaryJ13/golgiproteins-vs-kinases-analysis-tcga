# Figure 1: R equivalent of fig1_exome.py
#
# Optional alternative to the Python script; it produces the same figure.
# Run it from the working folder built by make_workdir.py, which already holds
# the input files:
#   glycosyltransferase.txt, Golgi_control_layer.txt, kinase.txt,
#   receptor_tyrosine_kinase.txt, exome_gene.csv, exome_class_dnds.csv,
#   cgc_by_class.csv
# From a terminal:   Rscript fig1_exome.R
# From RStudio:      setwd("path/to/work"); source("fig1_exome.R")
#
# Packages needed (once):
#   install.packages(c("ggplot2", "patchwork"))
#
# Writes Figure1.pdf and Figure1.png. Same file names as the Python script, so running both
# in the same folder overwrites the first output with the second.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- palette and style ---------------------------------------------------
ORD   <- c("receptor tyrosine kinase", "kinase", "Golgi control layer", "glycosyltransferase")
COL   <- c("receptor tyrosine kinase" = "#eda100", "kinase" = "#eb6834",
           "Golgi control layer" = "#1baf7a", "glycosyltransferase" = "#2a78d6")
SHAPE <- c("receptor tyrosine kinase" = 23, "kinase" = 24,
           "Golgi control layer" = 22, "glycosyltransferase" = 21)
LONG  <- c("receptor tyrosine kinase" = "Receptor tyrosine kinases",
           "kinase" = "Protein kinases",
           "Golgi control layer" = "Golgi control layer",
           "glycosyltransferase" = "Glycosylation enzymes")
SHORT <- c("receptor tyrosine kinase" = "RTKs", "kinase" = "Kinases",
           "Golgi control layer" = "Golgi control",
           "glycosyltransferase" = "Glycosylation")
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#dcdbd6"; RED <- "#c01c28"

# ggplot2 >= 3.5 needs legend.position.inside for legends placed inside the plot
legend_inside <- function(pos, just) {
  if (packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = pos,
          legend.justification = just)
  else
    theme(legend.position = pos, legend.justification = just)
}

base_theme <- theme_minimal(base_size = 8) +
  theme(plot.background  = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(colour = INK2, size = 8),
        axis.title       = element_text(colour = INK, size = 9),
        plot.title       = element_text(colour = INK, size = 10, face = "bold",
                                        hjust = 0, margin = ggplot2::margin(b = 6)),
        plot.title.position = "plot",
        legend.text      = element_text(colour = INK, size = 8),
        legend.title     = element_blank())

# ---- data --------------------------------------------------------------
# assign each gene to the first class it appears in (like setdefault in Python)
cls <- character(0)
for (c in ORD) {
  genes <- trimws(readLines(paste0(gsub(" ", "_", c), ".txt"), warn = FALSE))
  genes <- genes[genes != "" & !(genes %in% names(cls))]
  cls[genes] <- c
}

g <- read.csv("exome_gene.csv", stringsAsFactors = FALSE)
g$n <- rowSums(g[, c("n_syn", "n_mis", "n_non", "n_spl")])
g <- g[g$n > 0, ]
g$classe <- unname(cls[g$gene_name])
sub <- g[!is.na(g$classe), ]
sub$q <- pmax(sub$qallsubs_cv, 1e-17)
sub$w <- pmin(pmax(sub$wmis_cv, 0.05), 30)
sub$classe <- factor(sub$classe, levels = ORD)

cd <- read.csv("exome_class_dnds.csv", stringsAsFactors = FALSE)
cd <- cd[cd$stat == "wall", ]
rownames(cd) <- cd$classe
cg <- read.csv("cgc_by_class.csv", stringsAsFactors = FALSE)
rownames(cg) <- cg$gene_class

NSIG  <- sum(g$qallsubs_cv < 0.05)
NTEST <- nrow(g)

# vertical position of the classes in panels B, C and D (RTK at the top)
ypos <- setNames(rev(seq_along(ORD)), ORD)

# ---- Panel A --------------------------------------------------------
draw_order <- c("kinase", "receptor tyrosine kinase", "Golgi control layer", "glycosyltransferase")
subA <- sub[order(match(as.character(sub$classe), draw_order)), ]

# labels: name, horizontal offset (+ right, - left), vertical offset
lab <- data.frame(
  gene = c("BRAF", "STK11", "ATM", "MAP2K4", "MAP3K1", "EGFR", "FGFR2", "ERBB2", "AKT1", "ERBB3"),
  dx   = c(-1, 1, 1, 1, -1, 1, 1, 1, 1, -1),
  dy   = c(0.2, 0.2, -0.6, 0.2, 0.2, 0.2, -0.6, 0.2, -0.6, -0.6))
lab <- merge(lab, sub[, c("gene_name", "w", "q")], by.x = "gene", by.y = "gene_name")
lab$x <- lab$w * ifelse(lab$dx > 0, 1.08, 1 / 1.08)
lab$y <- -log10(lab$q) + lab$dy
lab$hjust <- ifelse(lab$dx > 0, 0, 1)

pA <- ggplot(subA, aes(w, -log10(q))) +
  geom_hline(yintercept = -log10(0.05), colour = RED, linewidth = 0.5, linetype = "22") +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.45, linetype = "32") +
  annotate("text", x = 0.055, y = -log10(0.05), label = " q = 0.05",
           colour = RED, size = 7.5 / .pt, hjust = 0, vjust = -0.4) +
  geom_point(aes(fill = classe, shape = classe), size = 2.1, colour = SURF,
             stroke = 0.3, alpha = 0.85) +
  geom_text(data = lab, aes(x, y, label = gene, hjust = hjust),
            size = 7 / .pt, colour = INK2, vjust = 0) +
  scale_fill_manual(values = COL, breaks = ORD, labels = LONG[ORD]) +
  scale_shape_manual(values = SHAPE, breaks = ORD, labels = LONG[ORD]) +
  scale_x_log10(limits = c(0.045, 40),
                breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10")) +
  scale_y_continuous(limits = c(-0.6, 26)) +
  labs(x = "dN/dS for missense substitutions", y = expression(-log[10]~q),
       title = "A    Selection acts on the kinases, not on the glycosylation machinery",
       subtitle = sprintf("TCGA MC3, 10,168 tumours; dNdScv fitted across %s genes exome-wide, %d significant",
                          format(NTEST, big.mark = ","), NSIG)) +
  base_theme +
  theme(plot.title = element_text(size = 10.5, face = "bold", colour = INK, hjust = 0),
        plot.subtitle = element_text(size = 8, colour = INK2, hjust = 0),
        axis.line = element_line(colour = GRID, linewidth = 0.3),
        legend.background = element_rect(fill = SURF, colour = GRID, linewidth = 0.3),
        legend.key.size = unit(10, "pt")) +
  legend_inside(c(0.01, 0.99), c(0, 1)) +
  guides(fill = guide_legend(override.aes = list(size = 2.6, colour = NA, alpha = 1)))

# ---- Panel B --------------------------------------------------------
bdat <- do.call(rbind, lapply(ORD, function(c) {
  s <- sub[sub$classe == c, ]; k <- sum(s$q < 0.05); n <- nrow(s)
  data.frame(classe = c, y = ypos[[c]], pct = 100 * k / n,
             lo = 100 * qbeta(0.025, k + 0.5, n - k + 0.5),   # intervallo di Jeffreys
             hi = 100 * qbeta(0.975, k + 0.5, n - k + 0.5),
             lab = sprintf("%d/%d", k, n))
}))

pB <- ggplot(bdat) +
  geom_rect(aes(xmin = 0, xmax = pct, ymin = y - 0.275, ymax = y + 0.275, fill = classe)) +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y), colour = INK, linewidth = 0.6) +
  geom_text(aes(x = hi + 0.9, y = y, label = lab), hjust = 0, size = 8 / .pt, colour = INK2) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(0, 27), expand = c(0, 0)) +
  scale_y_continuous(breaks = ypos, labels = SHORT[names(ypos)], limits = c(0.5, 4.7)) +
  labs(x = "Under positive selection,\nq < 0.05 (% of class)", y = NULL,
       title = "B    dN/dS, gene by gene") +
  base_theme +
  theme(panel.grid.major.y = element_blank(),
        axis.line.x = element_line(colour = GRID, linewidth = 0.3),
        axis.text.y = element_text(size = 8.5), axis.title.x = element_text(size = 8.5))

# ---- Panel C --------------------------------------------------------
cdat <- do.call(rbind, lapply(ORD, function(c) {
  r <- cg[c, ]; n <- r$n_class
  seq <- 100 * r$n_seq / n; oth <- 100 * (r$n_census - r$n_seq) / n
  data.frame(classe = c, y = ypos[[c]], seq = seq, tot = seq + oth,
             lab = sprintf("%d/%d", r$n_census, n))
}))

pC <- ggplot(cdat) +
  geom_rect(aes(xmin = seq, xmax = tot, ymin = y - 0.275, ymax = y + 0.275,
                fill = classe, alpha = "translocation or\ncopy number only")) +
  geom_rect(aes(xmin = 0, xmax = seq, ymin = y - 0.275, ymax = y + 0.275,
                fill = classe, alpha = "point mutation")) +
  geom_text(aes(x = tot + 1.4, y = y, label = lab), hjust = 0, size = 8 / .pt, colour = INK2) +
  scale_fill_manual(values = COL, guide = "none") +
  scale_alpha_manual(values = c("point mutation" = 1, "translocation or\ncopy number only" = 0.32),
                     breaks = c("point mutation", "translocation or\ncopy number only")) +
  guides(alpha = guide_legend(override.aes = list(fill = INK2))) +
  scale_x_continuous(limits = c(0, 62), expand = c(0, 0)) +
  scale_y_continuous(breaks = ypos, labels = NULL, limits = c(0.5, 4.7)) +
  labs(x = "In the Cancer Gene Census\n(% of class)", y = NULL,
       title = "C    COSMIC Census v104") +
  base_theme +
  theme(panel.grid.major.y = element_blank(),
        axis.line.x = element_line(colour = GRID, linewidth = 0.3),
        axis.title.x = element_text(size = 8.5),
        legend.text = element_text(size = 7, colour = INK),
        legend.background = element_rect(fill = SURF, colour = GRID, linewidth = 0.3),
        legend.key.size = unit(9, "pt")) +
  legend_inside(c(0.98, 0.03), c(1, 0))

# ---- Panel D --------------------------------------------------------
ddat <- data.frame(classe = ORD, y = ypos[ORD],
                   dNdS = cd[ORD, "dNdS"], lo = cd[ORD, "lo"], hi = cd[ORD, "hi"], p = cd[ORD, "p"])
ddat$lab <- paste0(sprintf("%.3f", ddat$dNdS), ifelse(ddat$p < 0.05, "", "   n.s."))

pD <- ggplot(ddat) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.5, linetype = "22") +
  annotate("text", x = 1, y = 4.55, label = "  neutral", hjust = 0, vjust = 0,
           size = 7 / .pt, colour = INK2) +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = classe),
               linewidth = 1.1, lineend = "round") +
  geom_point(aes(dNdS, y, fill = classe), shape = 21, size = 2.6, colour = SURF, stroke = 0.5) +
  geom_text(aes(x = hi + 0.003, y = y, label = lab), hjust = 0, size = 8 / .pt, colour = INK2) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(0.982, 1.195), breaks = seq(1, 1.175, 0.025),
                     labels = sprintf("%.3f", seq(1, 1.175, 0.025))) +
  scale_y_continuous(breaks = ypos, labels = SHORT[names(ypos)], limits = c(0.5, 4.8)) +
  labs(x = "dN/dS pooled over the whole class", y = NULL,
       title = "D    Pooled over the class, glycosylation is indistinguishable from neutral") +
  base_theme +
  theme(panel.grid.major.y = element_blank(),
        axis.line.x = element_line(colour = GRID, linewidth = 0.3),
        axis.text.y = element_text(size = 8.5))

# ---- assembly and saving -------------------------------------------
fig <- pA / (pB | pC) / pD +
  plot_layout(heights = c(1.42, 0.86, 0.74)) &
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        plot.margin = ggplot2::margin(6, 8, 6, 6))

fig

ggsave("Figure1.pdf", fig, width = 7.2, height = 10.2, units = "in", bg = SURF)
ggsave("Figure1.png", fig, width = 7.2, height = 10.2, units = "in", dpi = 300, bg = SURF)

cat("significativi esoma:", NSIG, "su", NTEST, "\n")
for (c in ORD) {
  s <- sub[sub$classe == c, ]
  cat(sprintf("  %-26s %d/%d\n", c, sum(s$q < 0.05), nrow(s)))
}
