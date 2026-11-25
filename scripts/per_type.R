# per_type.R -- dN/dS within each tumour type, exome-wide fit in each
#
# Answers the objection that pooling 33 tumour types hides selection confined
# to one of them. Each tumour type gets its own dNdScv run with
# gene_list = NULL, so the background is fitted exome-wide within that type,
# and the four gene classes are then estimated against that fit.
#
# Input : mut_tipo.tsv  -- run assign_tumour_type.py first; it labels every
#                          mutation in exome_filt.tsv with its TCGA tumour type
#         the four class .txt lists
# Output: per_type_genes.csv, per_type_class.csv
#
# Runtime ~90 min for the 11 types above 50,000 substitutions.

library(dndscv)

if (!file.exists("mut_tipo.tsv")) {
  stop("mut_tipo.tsv not found. Run first:  python assign_tumour_type.py")
}
m <- read.table("mut_tipo.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
classi <- c("glycosyltransferase", "Golgi_control_layer",
            "kinase", "receptor_tyrosine_kinase")
liste <- lapply(classi, function(c) readLines(paste0(c, ".txt")))
names(liste) <- classi

conteggi <- sort(table(m$tipo), decreasing = TRUE)
tipi <- names(conteggi)[conteggi >= 50000]
cat("tumour types to run:", length(tipi), "\n")

geni <- list(); classe <- list()

for (tp in tipi) {
  sub <- m[m$tipo == tp, c("sampleID", "chr", "pos", "ref", "mut")]
  cat("\n=====", tp, nrow(sub), "mutations,",
      length(unique(sub$sampleID)), "tumours =====\n")
  res <- try(dndscv(sub, gene_list = NULL,
                    max_muts_per_gene_per_sample = 3,
                    max_coding_muts_per_sample = Inf,
                    outp = 3, outmats = TRUE), silent = TRUE)
  if (inherits(res, "try-error")) { cat("FAILED:", res, "\n"); next }

  g <- res$sel_cv
  g$tipo <- tp
  geni[[tp]] <- g[, c("tipo", "gene_name", "n_syn", "n_mis", "n_non", "n_spl",
                      "wmis_cv", "wnon_cv", "pallsubs_cv", "qallsubs_cv")]

  for (cl in classi) {
    gl <- liste[[cl]][liste[[cl]] %in% g$gene_name]
    gs <- try(genesetdnds(res, gene_list = gl)$globaldnds_geneset, silent = TRUE)
    if (inherits(gs, "try-error")) next
    classe[[paste(tp, cl)]] <- data.frame(
      tipo = tp, classe = gsub("_", " ", cl), n_geni = length(gl),
      stat = sub("_geneset$", "", rownames(gs)),
      dNdS = signif(gs$mle, 4), lo = signif(gs$cilow, 4),
      hi = signif(gs$cihigh, 4), p = signif(gs$pval.lrt, 3))
  }
  cat("  significant genes (q<0.05):", sum(g$qallsubs_cv < 0.05, na.rm = TRUE), "\n")
  rm(res); gc()
}

write.csv(do.call(rbind, geni), "per_type_genes.csv", row.names = FALSE)
cc <- do.call(rbind, classe)
cc <- cc[cc$stat %in% c("wmis", "wall"), ]
write.csv(cc, "per_type_class.csv", row.names = FALSE)
cat("\ndone\n")
