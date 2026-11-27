# run_dndscv.R
#
# Exome-wide dN/dS on the TCGA MC3 call set, then per-class pooled estimates.
#
# Input : exome_filt.tsv  -- sampleID, chr, pos, ref, mut (with header)
#         (extract_all.py, then filter_hypermutators.py)
#         glycosyltransferase.txt, Golgi_control_layer.txt,
#         kinase.txt, receptor_tyrosine_kinase.txt -- one HGNC symbol per line
#
# Output: exome_gene.csv        -- dNdScv sel_cv, all 20,091 reference genes
#         exome_globaldnds.csv  -- pooled dN/dS over the whole exome
#         exome_class_dnds.csv  -- pooled dN/dS per gene class
#
# Requires: dndscv (Martincorena et al. 2017), with its bundled
#           RefCDS_human_GRCh37_dNdScv.0.1.0.rda reference.
#
# Runtime ~25 min, peak memory ~5 GB.

library(dndscv)

mut <- read.table("exome_filt.tsv", sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE,
                  col.names = c("sampleID", "chr", "pos", "ref", "mut"))
# expected: 2,096,036 substitutions in 10,168 tumours

cat("mutations:", nrow(mut), " samples:", length(unique(mut$sampleID)), "\n")

# --- the one decision that matters ---------------------------------------
# gene_list is deliberately left NULL. dNdScv fits its neutral background
# (negative binomial regression of synonymous counts on gene covariates)
# from whatever genes it is given. Restricting it to the 964 genes of
# interest makes it learn a baseline from a functionally related,
# non-random set, which produced significant *negative* selection on
# three genes (CHST4, CHST6, GOPC). Fitting across all genes and 
# applying the class labels afterwards removes it.
dnds <- dndscv(mut,
               gene_list = NULL,
               max_muts_per_gene_per_sample = 3,
               max_coding_muts_per_sample = Inf,  # already filtered upstream
               outp = 3,
               outmats = TRUE)   # required by genesetdnds below

write.csv(dnds$sel_cv,     "exome_gene.csv",       row.names = FALSE)
write.csv(dnds$globaldnds, "exome_globaldnds.csv", row.names = FALSE)

# The run above costs about 50 minutes. Save the object so that a mistake in
# the class step below does not mean paying for it twice.
saveRDS(dnds, "dndsout.rds", compress = FALSE)

# --- pooled dN/dS per class ----------------------------------------------
# genesetdnds re-estimates the ratio over a gene set using the background
# model already fitted above, so the classes share one common baseline.
classes <- c("glycosyltransferase", "Golgi_control_layer",
             "kinase", "receptor_tyrosine_kinase")

out <- do.call(rbind, lapply(classes, function(cl) {
  genes <- readLines(paste0(cl, ".txt"))
  genes <- genes[genes %in% dnds$sel_cv$gene_name]
  # NB: the element is `globaldnds_geneset`, not `globaldnds`, and the
  # parameter name lives in the row names, not in a column. Its p value
  # column is `pval.lrt` -- the likelihood ratio test, which the package
  # recommends over the Wald p value also returned.
  g <- genesetdnds(dnds, gene_list = genes)$globaldnds_geneset
  data.frame(classe = gsub("_", " ", cl),
             stat   = sub("_geneset$", "", rownames(g)),
             dNdS   = signif(g$mle, 4),
             lo     = signif(g$cilow, 4),
             hi     = signif(g$cihigh, 4),
             p      = signif(g$pval.lrt, 3))
}))

out <- out[out$stat %in% c("wmis", "wnon", "wall"), ]
write.csv(out, "exome_class_dnds.csv", row.names = FALSE)

cat("\ndone\n")
print(out)
