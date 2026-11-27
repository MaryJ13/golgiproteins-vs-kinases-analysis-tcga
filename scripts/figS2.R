# Supplementary Figure S2: R equivalent of figS2.py
#
# Optional alternative to the Python script; it produces the same figure.
# Run it from the working folder built by make_workdir.py, which already holds
# the input files:
#   per_type_class.csv, genemuts.csv, exome_gene.csv,
#   glycosyltransferase.txt, Golgi_control_layer.txt, kinase.txt,
#   receptor_tyrosine_kinase.txt
# From a terminal:   Rscript figS2.R
# From RStudio:      setwd("path/to/work"); source("figS2.R")
#
# Packages needed (once):
#   install.packages(c("ggplot2", "patchwork"))
#
# Writes FigureS2.pdf and FigureS2.png, the same file names as the Python script,
# so running both in one folder overwrites the first output with the second.

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
AB <- c("Skin Cutaneous Melanoma" = "SKCM", "Lung adenocarcinoma" = "LUAD",
        "Lung squamous cell carcinoma" = "LUSC",
        "Uterine Corpus Endometrial Carcinoma" = "UCEC",
        "Colon adenocarcinoma" = "COAD", "Stomach adenocarcinoma" = "STAD",
        "Bladder Urothelial Carcinoma" = "BLCA", "Breast invasive carcinoma" = "BRCA",
        "Head and Neck squamous cell carcinoma" = "HNSC",
        "Cervical squamous cell carcinoma and endocervical adenocarcinoma" = "CESC",
        "Liver hepatocellular carcinoma" = "LIHC")
NOTE <- "#555555"

# ggplot2 >= 3.5 needs legend.position.inside for legends placed inside the plot
legend_inside <- function(pos, just) {
  if (packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = pos,
          legend.justification = just)
  else
    theme(legend.position = pos, legend.justification = just)
}

base_theme <- theme_classic(base_size = 8) +
  theme(axis.line   = element_line(colour = "#444444", linewidth = 0.4),
        axis.ticks  = element_line(colour = "#444444", linewidth = 0.4),
        axis.text   = element_text(colour = "black", size = 7.5),
        axis.title  = element_text(size = 8),
        plot.title  = element_text(size = 9, face = "bold", hjust = 0),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.text  = element_text(size = 6.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.caption = element_text(size = 6.6, colour = NOTE, hjust = 0),
        plot.caption.position = "panel")

# ---- Panel A --------------------------------------------------------
cc <- read.csv("per_type_class.csv", stringsAsFactors = FALSE)
w  <- cc[cc$stat == "wall", ]

# tumour types ordered by kinase dN/dS, highest at the top
wk   <- w[w$classe == "kinase", ]
tipi <- wk$tipo[order(-wk$dNdS)]
OFF  <- setNames(c(0.27, 0.09, -0.09, -0.27), ORD)

w <- w[w$tipo %in% tipi & w$classe %in% ORD, ]
w$i      <- match(w$tipo, tipi) - 1
w$y      <- (length(tipi) - 1 - w$i) - OFF[w$classe]   # same order as the inverted axis in Python
w$sig    <- ifelse(w$p < 0.05, "p < 0.05", "n.s.")
w$stroke <- ifelse(w$p < 0.05, 0.55, 0)
w$size   <- ifelse(w$p < 0.05, 1.9, 1.55)
w$classe <- factor(w$classe, levels = ORD)

pA <- ggplot(w) +
  geom_vline(xintercept = 1, colour = "#888888", linewidth = 0.45, linetype = "dashed")
for (cl in ORD) {   # one level per class, so the bar colour stays fixed
  pA <- pA + geom_segment(data = w[w$classe == cl, ],
                          aes(x = lo, xend = hi, y = y, yend = y),
                          colour = COL[[cl]], linewidth = 0.6, alpha = 0.85)
}
pA <- pA +
  geom_point(aes(dNdS, y, fill = classe, shape = classe, colour = sig,
                 size = size, stroke = stroke)) +
  scale_size_identity() +
  scale_fill_manual(values = COL, breaks = ORD, labels = LONG[ORD]) +
  scale_shape_manual(values = SHAPE, breaks = ORD, labels = LONG[ORD]) +
  scale_colour_manual(values = c("p < 0.05" = "#222222", "n.s." = "#222222"),
                      breaks = "p < 0.05") +
  guides(fill   = guide_legend(order = 1, nrow = 1,
                               override.aes = list(size = 2, stroke = 0)),
         shape  = guide_legend(order = 1, nrow = 1),
         colour = guide_legend(order = 2,
                               override.aes = list(shape = 21, fill = "#cccccc",
                                                   size = 2, stroke = 0.6))) +
  scale_x_continuous(limits = c(0.64, 1.78), expand = c(0, 0),
                     breaks = seq(0.8, 1.6, 0.2)) +
  scale_y_continuous(breaks = rev(seq_along(tipi)) - 1,
                     labels = ifelse(tipi %in% names(AB), AB[tipi], tipi),
                     limits = c(-0.6, length(tipi) - 0.4), expand = c(0, 0)) +
  labs(x = "dN/dS pooled over the class, within each tumour type", y = NULL,
       title = "A   Selection is confined to the kinases in every tumour type",
       caption = paste0(
         "11 tumour types with >50,000 substitutions; each fitted exome-wide in its own run. No Golgi or\n",
         "glycosylation gene reaches significance in any single type (0/11 for both classes, and 0/11 at gene level).")) +
  base_theme +
  theme(panel.grid.major.x = element_line(colour = "#eeeeee", linewidth = 0.35),
        legend.position = "top", legend.justification = "center",
        legend.box = "horizontal", legend.spacing.x = grid::unit(2, "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.key.width = grid::unit(10, "pt"))

# ---- Panel B --------------------------------------------------------
gm  <- read.csv("genemuts.csv", stringsAsFactors = FALSE)
sel <- read.csv("exome_gene.csv", stringsAsFactors = FALSE)
d <- merge(gm, sel[, c("gene_name", "wmis_cv", "pmis_cv", "qallsubs_cv")], by = "gene_name")
d <- d[d$exp_mis > 0, ]
pm <- d$pmis_cv; pm[is.na(pm)] <- 1
d$qmis <- p.adjust(pm, method = "BH")        # = multipletests(..., "fdr_bh")
det <- d[d$qmis < 0.05 & d$wmis_cv > 1, ]
det <- det[order(det$exp_mis), ]

# empirical detection boundary: minimum dN/dS detected in a moving window
xs <- 10^seq(log10(min(det$exp_mis)), log10(max(det$exp_mis)), length.out = 60)
ys <- sapply(xs, function(x) {
  v <- det$wmis_cv[det$exp_mis >= x / 1.6 & det$exp_mis <= x * 1.6]
  if (length(v)) min(v) else NA
})
fill_na <- function(v) {                     # forward-fill then back-fill, as in pandas
  for (k in seq_along(v)[-1]) if (is.na(v[k])) v[k] <- v[k - 1]
  for (k in rev(seq_along(v))[-1]) if (is.na(v[k])) v[k] <- v[k + 1]
  v
}
ys <- fill_na(ys)
ys <- sapply(seq_along(ys), function(k)      # centred moving average, window of 7
  mean(ys[max(1, k - 3):min(length(ys), k + 3)]))
bound <- data.frame(x = xs, y = ys)

# classes: if a gene is in several files, the last one read wins (like the Python dict)
classi <- character(0)
for (cfile in c("glycosyltransferase", "Golgi_control_layer", "kinase",
                "receptor_tyrosine_kinase")) {
  genes <- scan(paste0(cfile, ".txt"), what = "", quiet = TRUE)
  classi[genes] <- gsub("_", " ", cfile)
}
d$classe <- unname(classi[d$gene_name])
d$ntot <- rowSums(d[, c("n_syn", "n_mis", "n_non", "n_spl")])
cl <- d[!is.na(d$classe) & d$ntot > 0, ]
cl$wplot <- pmin(pmax(cl$wmis_cv, 0.3), 30)
cl$classe <- factor(cl$classe, levels = ORD)
cl <- cl[order(cl$classe), ]                  # draw in the order of ORD

labB <- cl[cl$gene_name %in% c("BRAF", "EGFR", "FGFR2", "B3GNT1", "YIPF6"), ]
labB$ylab <- pmin(labB$wmis_cv, 30)

noteB <- grid::textGrob(
  paste0("Red line: the weakest missense dN/dS that genes of comparable expectation actually reached significance with.\n",
         "Shaded: the detectable region. No glycosylation enzyme enters it; the highest is B3GNT1 at 1.95 (q = 0.43)."),
  x = grid::unit(0.012, "npc"), y = grid::unit(0.022, "npc"), hjust = 0, vjust = 0,
  gp = grid::gpar(fontsize = 6.6, col = NOTE))

pB <- ggplot() +
  geom_ribbon(data = bound, aes(x = x, ymin = y, ymax = 40), fill = "#f3b0a0", alpha = 0.30) +
  geom_hline(yintercept = 1, colour = "#888888", linewidth = 0.45, linetype = "dashed") +
  geom_point(data = cl, aes(exp_mis, wplot, fill = classe, shape = classe),
             size = 1.15, stroke = 0, alpha = 0.62) +
  geom_line(data = bound, aes(x, y), colour = "#c0392b", linewidth = 0.75) +
  geom_text(data = labB, aes(exp_mis * 1.07, ylab * 1.06, label = gene_name),
            hjust = 0, vjust = 0, size = 6.6 / .pt, colour = "#333333") +
  annotation_custom(noteB) +
  scale_fill_manual(values = COL, breaks = ORD, labels = LONG[ORD]) +
  scale_shape_manual(values = SHAPE, breaks = ORD, labels = LONG[ORD]) +
  guides(fill = guide_legend(override.aes = list(size = 1.8, alpha = 0.8))) +
  scale_x_log10(breaks = c(10, 100, 1000),
                labels = c(expression(10^1), expression(10^2), expression(10^3))) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20), labels = c("0.5", "1", "2", "5", "10", "20"),
                limits = c(0.28, 48), expand = c(0, 0)) +
  annotation_logticks(sides = "bl", colour = "#444444",
                      short = grid::unit(0.6, "mm"), mid = grid::unit(0.9, "mm"), long = grid::unit(1.2, "mm")) +
  labs(x = "Missense substitutions expected under neutrality (gene size x burden)",
       y = "dN/dS for missense substitutions",
       title = "B   What the analysis could have detected",
  caption = paste0(
    "Red line: the weakest missense dN/dS that genes of comparable expectation actually reached significance with.\n",
    "Shaded: the detectable region. No glycosylation enzyme enters it; the highest is B3GNT1 at 1.95 (q = 0.43).")) +
  base_theme +
  theme(panel.grid.major = element_line(colour = "#eeeeee", linewidth = 0.35),
        legend.background = element_rect(fill = "#FFFFFFF2", colour = "#dddddd",
                                         linewidth = 0.3),
        legend.key.size = grid::unit(9, "pt")) +
  legend_inside(c(0.99, 0.99), c(1, 1))

# ---- assembly and saving -----------------------------------------------
fig <- pA / pB + plot_layout(heights = c(1.15, 1)) &
  theme(plot.margin = ggplot2::margin(6, 10, 6, 6))

fig

ggsave("FigureS2.pdf", fig, width = 7.2, height = 8.4, units = "in")
ggsave("FigureS2.png", fig, width = 7.2, height = 8.4, units = "in", dpi = 220)
cat("written FigureS2.pdf / .png\n")
