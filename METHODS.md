# Methods — Figure 1 and Supplementary Figure S1

Selection analysis of the glycosylation machinery, the Golgi control layer,
protein kinases and receptor tyrosine kinases across TCGA.

What was done, in order, with the parameters and the reasons.

---

## 1. Rationale for the measure

The question is whether the glycosylation machinery shows evidence of
positive selection in cancer genomes, as protein kinases do.

The direct approach, comparing the proportion of tumours carrying a mutation
in each class, cannot answer it. Mutation frequency scales with coding
length, and the classes differ systematically in length: the median protein
kinase has 1,860 coding bases and the longest, TTN, has 108,000, while
glycosyltransferases have a median of 1,332 and none above 3,800. Run that
way, glycosylation enzymes appeared four times less likely than random genes
to be mutated in ≥2% of tumours (p = 0.004). Refitting with coding length as
an offset in a negative binomial model removed the effect entirely:
length-adjusted rate ratios were 1.02 for glycosylation enzymes and 0.99 for
kinases, neither different from random. The frequency comparison measures
gene size.

dN/dS has no such confounding. It is the ratio of observed to expected
protein-altering substitutions, where the expectation is set by the
synonymous substitutions **in the same gene**. Coding length, trinucleotide
composition and the tumour's overall mutational burden enter numerator and
denominator identically and cancel. A value of 1 indicates neutrality, above
1 positive selection.

## 2. Mutation data

Somatic mutations were taken from the TCGA MC3 public call set
(`mc3.v0.2.8.PUBLIC.maf.gz`, NCI Genomic Data Commons, GRCh37), the
consortium-harmonised pan-cancer variant set across 33 tumour types.

Processing (`extract_all.py`):

1. Single-base substitutions only; reference and alternate alleles both in
   {A, C, G, T} and different from each other. Indels and dinucleotide or
   complex events were discarded.
2. Chromosome identifiers normalised. **MC3 encodes X as `23` and Y as `24`**,
   which dNdScv does not recognise; left uncorrected, 5,937 mutations are
   silently dropped, including those in *ATP7A*, *SLC35A2*, *C1GALT1C1* and
   *PIGA*.
3. Sample barcodes truncated to the 15-character tumour identifier.
4. Exact duplicates (same sample, position, ref and alt) removed.

This yielded 3,284,873 substitutions in 10,294 tumours.

Hypermutated tumours were then removed (`filter_hypermutators.py`): 126
samples carrying more than 3,000 coding substitutions, the MSI and POLE-mutant
cases. dNdScv fits one background mutation model across all samples, and
samples carrying a thousand times the median distort it; they also dominate
per-gene counts. The threshold is not delicate — the distribution is strongly
bimodal and any cut between 2,000 and 5,000 removes the same cluster. The
analysis set is **2,096,036 substitutions in 10,168 tumours**.

## 3. Gene classes

Four mutually exclusive classes.

**Glycosylation enzymes and the Golgi control layer** were taken from a
curated in-house gene set (GolgiSet, 507 genes). Rows annotated
`Glycosylation enzyme` form the first class; membrane trafficking components,
channels and transporters, and lipid trafficking and metabolism genes form
the second.

**Protein kinases and receptor tyrosine kinases** were taken from KinHub, the
Manning kinome with current gene symbols. Receptor tyrosine kinases are the
TK group excluding the families whose members are cytoplasmic (Abl, Ack, Csk,
FAK, Fer, Jak, Src, Syk, Tec, TK-Unique) and two receptor families that are
pseudokinases (CCK4/PTK7, Lmr/LMTK); *ERBB3* is a receptor of the EGFR family
but is catalytically dead and was assigned to the kinase class. Before
filtering this rule returns exactly 58 receptor tyrosine kinases, the
canonical count across their 20 families, which serves as the check that the
KinHub table was parsed correctly.

Three filters were then applied, in order:

1. **Absent from the dNdScv GRCh37 reference** — 6 GolgiSet genes and 14
   kinome entries. A gene the reference does not know has no dN/dS. Where a
   gene exists in the reference under a superseded symbol, an explicit alias
   table maps it (e.g. `B3GLCT`→`B3GALTL`, `GRK2`→`ADRBK1`); every
   substitution is logged.
2. **Present in both the GolgiSet and the kinome** — *CERT1*, *POMK*,
   *SCYL1*, *SCYL3*. These were retained in the Golgi classes and removed
   from the kinase classes, so that no gene belonging to the object under
   study also serves as its control.
3. **Not a catalytically active protein kinase** — 26 entries in KinHub's
   Atypical and Other groups: ubiquinone kinases (*ADCK1/2/5*),
   pseudokinases with degenerate catalytic residues (*TRIB1/2/3*,
   *STRADA/B*, *MLKL*, *PXK*, *RYK*, *AATK*, *SCYL2*, *ULK4*, *RNASEL*,
   *RPS6KC1*, *RPS6KL1*, *IRAK2*), and scaffold or chromatin proteins
   carrying a kinase-like domain (*BAZ1A*, *TRIM33*, *TRIM66*, *TRRAP*,
   *PAN3*, *ABR*), plus *STK19* and *PRPF4B*. This matters because the
   kinase classes are the positive control: several of these genes are
   COSMIC Census cancer genes, and leaving them in would inflate the control
   with signal unrelated to kinase activity. The filter makes the control
   harder to pass, not easier.

Finally, genes with no observed coding substitution were excluded from the
figure, since dNdScv returns dN/dS = 0 for them, which is not an estimate.
Three were affected: *OST4*, *GOLGA7B* and *SBK3*.

Final class sizes: **230** glycosylation enzymes, **269** Golgi control layer
genes, **434** protein kinases, **51** receptor tyrosine kinases.

## 4. Selection analysis

dN/dS was estimated with dNdScv 0.1.0 (Martincorena et al., *Cell* 2017)
against its bundled GRCh37 coding reference, under R 4.3.3. dNdScv models
substitution rates with a 192-parameter trinucleotide context model and
regresses synonymous counts on gene-level covariates (expression,
replication timing, chromatin state) by negative binomial regression, which
supplies each gene's expected neutral rate.

**The analysis was run once, across all 20,091 reference genes
(`gene_list = NULL`), with the class labels applied afterwards.** This is the
single methodological decision that most affects the result, and it was
reached only after two other designs produced artefacts:

| Design | Background fitted on | Glyco | Golgi | Kinases | RTKs |
|---|---|---|---|---|---|
| One run per class | each class alone | 0/225 | 0/264 | 34/424 | 8/51 |
| One run, class genes only | those 964 genes | 4/225 | 2/264 | 31/424 | 7/51 |
| **One run, whole exome** | **19,155 genes** | **0/230** | **0/269** | **13/434** | **6/51** |

Under the second design three sulfotransferases reached significance with
dN/dS *below* 1 — *CHST4* at 0.55 (q = 0.015), *CHST6* at 0.79 (q = 0.050),
*GOPC* at 0.85 (q = 0.010) — that is, apparent significant negative selection
across ten thousand tumours. That cannot be real, and it is the signature of
a background model learned from a functionally related, non-random set of
genes. Under the exome-wide fit all three fall away (q = 0.21, 0.36, 0.23).

Multiple testing was controlled by Benjamini–Hochberg across all 19,155
tested genes (`qallsubs_cv`), not within each class. This is the conservative
choice and the reason *KIT* and *PDGFRA*, both established drivers, do not
reach significance here: the threshold the glycosylation enzymes fail to
clear is the same one *KIT* fails to clear.

Class-level pooled estimates (panel D) were obtained with `genesetdnds`
against that **same** exome-wide fit — one call per class — so that all four
panels share one background. Note that passing `gene_list` to `dndscv()`
instead would subset the coding reference before annotation and estimate each
class's substitution model from its own mutations alone; doing so changes the
Golgi control layer from 1.023 (p = 0.154) to 1.043 (p = 0.014), i.e. it
manufactures significant positive selection where there is none. `genesetdnds`
requires `outmats = TRUE` on the parent run.

## 5. Cancer Gene Census

Census membership was taken from COSMIC release v104 (19 May 2026), matched
on gene symbol. The `MUTATION_TYPES` field was decoded as: point mutation
(Mis, N, F, S, O), translocation (T), amplification (A), large deletion (D).
A gene may carry several, so the mechanism counts are not mutually exclusive.

The code letters were verified against genes whose mechanism is known
independently rather than taken on trust: *MDM2* carries A alone, *MYC* A and
T, *EGFR* A, O and Mis, *TP53* Mis, N, F, T and D.

Amplification and large deletion are reported in separate columns. They are
both copy-number events but point in opposite directions — a dosage-driven
oncogene against a deleted tumour suppressor — and among the kinases the two
are equally common (6 genes each), so a combined "copy number" column is
misleading.

## 6. Figure construction

**Panel A**, volcano: `wmis_cv` against −log10 `qallsubs_cv`, one point per
gene, 984 class genes. Values are clipped to [0.05, 30] on the x axis and
floored at 1e-17 on the y axis for display only; clipping never affects
whether a point counts as significant.

**Panel B**, the fraction of each class reaching q < 0.05, with 95% Jeffreys
binomial intervals, `Beta(k+0.5, n−k+0.5)` quantiles. Jeffreys rather than
Wald because two classes have k = 0, where a Wald interval collapses to zero
width and conveys nothing.

**Panel C**, Census membership as a percentage of class, split into genes
carrying at least one point-mutation code and the remainder. No error bars:
this is a count from a catalogue, not an estimate.

**Panel D**, pooled dN/dS per class (`wall` of `globaldnds_geneset`) with 95%
confidence intervals, marked n.s. where the likelihood ratio p exceeds 0.05.

**Supplementary Figure S1** carries the two checks of Sections 8 and 9:
(A) pooled dN/dS per class within each tumour type, with 95% confidence
intervals and outlined markers where p < 0.05; (B) observed missense dN/dS
against expected neutral missense count, one point per class gene, with the
empirical detection boundary drawn as a line and the detectable region
shaded.

## 7. Validation of the method

The analysis was run blind to the gene classes across the genome. Ranked by
p value, its twelve most significant genes are *TP53*, *PTEN*, *FBXW7*,
*KEAP1*, *NRAS*, *PIK3CA*, *NFE2L2*, *BRAF*, *KRAS*, *CTNNB1*, *IDH1* and
*HRAS* — the canonical pan-cancer driver set, with the expected effect sizes
(*KRAS* dN/dS 26.5, *NRAS* 24.4, *TP53* 15.2, *PIK3CA* 12.2, *BRAF* 10.5).
Thirty-six genes reach q = 0 to machine precision, so the ordering among the
very top is arbitrary and should not be read as a ranking. 145 genes of the
19,155 carrying at least one substitution reach q < 0.05, the usual yield for
this analysis on TCGA.

Independent re-execution of the full pipeline from the downloaded MAF
reproduced `exome_gene.csv` byte for byte and the global dN/dS estimates to
the last reported digit.

## 8. Detection power

A negative result is interpretable only against the effect size the test could
have found, so the detection threshold was calibrated on the analysis itself
rather than assumed.

Among all tested genes, 61 reached missense-driven significance with dN/dS
above 1 (Benjamini–Hochberg on dNdScv's missense p value; genes below 1 were
excluded so that tumour suppressors detected through truncating mutations do
not set the bar). For each gene, the detection threshold is the weakest
missense dN/dS reached by any of those 61 whose expected neutral missense
count lies within a factor of 1.3 of its own.

The calibration is checked on the receptor tyrosine kinases, where it predicts
exactly the six genes that are in fact significant.

A Poisson power calculation was tried first and discarded. It treats the
neutral expectation as known, whereas dNdScv fits it by negative binomial
regression whose dispersion adds real uncertainty, and it was visibly
anti-conservative: five non-significant genes sat above its threshold. The
empirical calibration has no such failure.

Median thresholds: 1.61 for receptor tyrosine kinases, 1.93 for protein
kinases, 1.98 for the Golgi control layer, 2.22 for glycosylation enzymes.
Per-gene values are in `results/power_by_gene.csv`.

## 9. Tumour-type stratification

Pooling 33 tumour types would dilute selection confined to one of them, so
each of the largest types was also analysed alone.

Tumours were assigned to a TCGA study from the tissue source site code, the
sixth and seventh characters of the barcode, against the GDC code table
(`source_data/tss.tsv`). All 10,168 tumours mapped, with no unrecognised
codes.

The 11 types carrying more than 50,000 substitutions, together holding 82.7%
of the analysis set, were each analysed in a separate dNdScv run with
`gene_list = NULL`, so that the background was fitted exome-wide *within* that
type; the four classes were then estimated against that fit with
`genesetdnds`. Per-type yields were 20 to 47 significant genes, consistent
with the reduced sample size.

This analysis is less powered than the pooled one, by construction. Its
purpose is not to find weaker signals but to catch a signal that pooling would
have hidden, which it would show as a strong effect in one type.

---

## References

Martincorena I, Raine KM, Gerstung M, Dawson KJ, Haase K, Van Loo P,
Davies H, Stratton MR, Campbell PJ. Universal patterns of selection in cancer
and somatic tissues. *Cell* 2017;171:1029–1041.

Ellrott K, Bailey MH, Saksena G, et al. Scalable open science approach for
mutation calling of tumor exomes using multiple genomic pipelines.
*Cell Systems* 2018;6:271–281. (The MC3 call set.)

Sondka Z, Dhir NB, Carvalho-Silva D, et al. COSMIC: a curated database of
somatic variants and clinical data for cancer. *Nucleic Acids Research*
2024;52:D1210–D1217. (The Cancer Gene Census.)

Manning G, Whyte DB, Martinez R, Hunter T, Sudarsanam S. The protein kinase
complement of the human genome. *Science* 2002;298:1912–1934. (KinHub.)
