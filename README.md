# Is the glycosylation machinery under positive selection in cancer?

Code and data for Figure 1 and Supplementary Figure S1 of:

> M. L. Maresca, R. Rizzo, F. Russo, S. Parashuraman and D. Russo
> **Tuning, not switching: the Golgi control layer as an analogue determinant of chemoresistance.**
> *Bioscience Reports*, invited review (submitted 2026).

Archived version: 

Repository authors: Maria Luigia Maresca (CNR, Naples; ORCID [0009-0006-4170-2422](https://orcid.org/0009-0006-4170-2422)), Riccardo Rizzo (CNR, Lecce; ORCID [0000-0002-5716-9195](https://orcid.org/0000-0002-5716-9195)), Francesco Russo (CNR, Naples; ORCID [0000-0002-7017-8344](https://orcid.org/0000-0002-7017-8344)), Seetharaman Parashuraman (CNR, Naples; ORCID [0000-0001-5113-4592](https://orcid.org/0000-0001-5113-4592)), Domenico Russo (CNR, Naples; ORCID [0000-0003-2171-657X](https://orcid.org/0000-0003-2171-657X))

---

## What we set out to test

The review argues that the glycosylation machinery contributes to cancer
through dosage and organisation rather than through mutation. That is a claim
about sequence: if these enzymes acted like oncogenes, selection would have
left a mark on their coding sequence, as it has on the kinases.

**Across ten thousand tumours, is there any signal of positive selection on the 
genes of the glycosylation machinery or of the Golgi control layer — searched 
with the method that finds *TP53* and *KRAS*?**

The analysis here shows any signal of positive selection, together with
the two checks that make a negative result interpretable: what effect size
the test could have detected, and whether pooling tumour types hid anything.

## What was done

Four gene classes — glycosylation enzymes, Golgi control layer, protein
kinases, receptor tyrosine kinases — tested for selection with dN/dS across
2,096,036 somatic substitutions from 10,168 TCGA tumours.

dN/dS is used rather than mutation frequency because frequency scales with
coding length, and these classes differ systematically in length. The naive
comparison appears to support the claim and does so for the wrong reason: the
apparent effect disappears entirely once length is accounted for.

**Headline result.** Of 145 genes under positive selection genome-wide, 6 are
receptor tyrosine kinases and 13 are protein kinases. None is a glycosylation
enzyme and none is a Golgi gene (0/230 and 0/269). The analysis would have
detected a 2.2-fold missense excess in a glycosylation enzyme of typical size;
the highest observed is 1.95. Analysed within each of the 11 largest tumour
types separately, kinases reach significance in 10 of 11 and the Golgi classes
in none.

**What the analysis cannot see, and why that matters.** dN/dS is blind to gene
dosage. A gene amplified in a fifth of tumours and never point-mutated returns
a ratio of exactly 1. *GOLPH3* is the worked example: unremarkable here
(dN/dS 1.18, q = 0.885) and nonetheless a COSMIC Census cancer gene whose only
recorded mutation type is amplification. Effects that act through gene dosage are 
therefore outside what this analysis tests.

`METHODS.md` gives the full method, with the reasoning behind each choice.

## Contents

```
README.md          this file
METHODS.md         the written method, in full
LICENSE.md         MIT for the code, CC BY 4.0 for data and text,
                   and the provenance of third-party material
CITATION.cff       how to cite this deposit
requirements.txt   Python packages (R packages: scripts/ENVIRONMENT.md)
make_workdir.py    run this first: builds the flat working folder
scripts/           the pipeline, in the order it runs, plus R versions
                   of the two figure scripts
source_data/       the inputs; two more are linked, not shipped
results/           everything the pipeline produced, and the two figures
supplementary/     Supplementary Tables S1-S7 of the review, as one workbook
```

### `source_data/`

Four of the six inputs are here. `GolgiSet.xlsx` is curated in-house and
exists nowhere else. `KinHub.pdf` is the published kinome as we parsed it,
kept as a file so the parse can be checked against the same table rather than
a later revision of the site. `tss.tsv` maps TCGA barcode tissue-source-site
codes to tumour types. And `exome_filt.tsv.gz` is the mutation table already
extracted and filtered: 2,096,036 substitutions in 10,168 tumours, 14 MB
gzipped, which saves a 700 MB download and about 35 minutes of processing.

Two inputs are linked rather than included:

| | Where | Why not here | Needed for |
|---|---|---|---|
| `mc3.v0.2.8.PUBLIC.maf.gz` | <https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc> | 700 MB, already public, versioned and permanent at the GDC | rebuilding `exome_filt.tsv`, which is already here |
| COSMIC Cancer Gene Census v104 | <https://cancer.sanger.ac.uk/cosmic> | its academic licence forbids redistribution | regenerating the Census annotation, which is already in `results/` |

**Neither is needed to reproduce any panel of either figure.** `get_mc3.py`
downloads the first automatically; the second has to be fetched by hand
because of the login. `source_data/README.md` gives the hashes, the exact
release to take, and which files circulating from earlier drafts are *not*
inputs.

### `results/`

Everything the pipeline produced, so that any number in the paper can be
checked without rerunning anything.

| File | What it holds |
|---|---|
| `Figure1_source_data.csv` | **open this first** — the 984 class genes, with counts, dN/dS, q and Census annotation |
| `exome_gene.csv` | dNdScv output for all 20,091 reference genes, unfiltered |
| `exome_class_dnds.csv` | pooled dN/dS per class, pan-cancer (Figure 1D) |
| `exome_globaldnds.csv` | pooled dN/dS over the whole exome |
| `genemuts.csv` | per-gene expected counts under neutrality; the input to the power calibration |
| `cgc_by_class.csv` | Cancer Gene Census counts by class and mechanism (Figure 1C) |
| `power_by_gene.csv` | per-gene detection threshold against observed dN/dS (Figure S1B) |
| `per_type_class.csv` | pooled dN/dS per class within each tumour type (Figure S1A) |
| `per_type_genes.csv.gz` | the full per-gene dNdScv output for each of the 11 tumour types; the evidence behind "no Golgi or glycosylation gene is significant in any single type" |
| `glycosyltransferase.txt`, `Golgi_control_layer.txt`, `kinase.txt`, `receptor_tyrosine_kinase.txt` | the four gene classes as used |
| `golgi_lists_log.txt`, `kinase_lists_log.txt` | every symbol substitution and every gene dropped, with the reason |
| `hypermutators.txt` | the 126 MSI and POLE tumours removed, with their substitution counts |
| `Figure1.pdf`, `FigureS1.pdf` | the figures as the scripts produced them |

### `supplementary/`

`Supplementary_Tables_S1-S7.xlsx` is the workbook of supplementary tables
submitted with the review (per-gene results, pooled dN/dS by class, Census
membership, background-model comparison, detection power, per-tumour-type
results). No script in this deposit regenerates the workbook itself.

## Three ways to run this

First build the working folder. The scripts read and write the current
directory, so they need everything in one place:

```bash
python make_workdir.py     # copies the three folders into ./work
cd work
python check_inputs.py     # verifies checksums, says what you can run
```

`check_inputs.py` prints which of the three routes below are open to you.

**1. Redraw the figures.** No downloads, seconds.

```bash
python fig1_exome.py      # -> Figure1.pdf
python figS1.py           # -> FigureS1.pdf
```

The same two figures can be drawn in R instead (needs `ggplot2` and
`patchwork`; both scripts write the same file names, so use one language or
the other in a given working folder):

```bash
Rscript fig1_exome.R      # -> Figure1.pdf, Figure1.png
Rscript figS1.R           # -> FigureS1.pdf, FigureS1.png
```

**2. Reanalyse from the filtered mutation table.** Still no downloads: this is
the full dN/dS analysis, starting from the 2,096,036 substitutions included
in `source_data/`.

```bash
Rscript run_dndscv.R           # ~50 min, ~6 GB RAM
python assign_tumour_type.py   # -> mut_tipo.tsv, needed by per_type.R
Rscript per_type.R             # ~90 min
python power2.py
python fig1_exome.py ; python figS1.py
```

**3. Rebuild everything from the original sources.** Adds the 700 MB download
and the extraction. The Census step additionally needs the COSMIC file; omit
it and the shipped annotation is used instead.

```bash
python get_mc3.py                                      # ~700 MB download
python extract_all.py
python filter_hypermutators.py mutations_all.tsv.gz
python build_kinase_lists.py KinHub.pdf
python build_golgi_lists.py GolgiSet.xlsx exome_gene.csv
python build_kinase_classes.py kinhub_parsed.csv exome_gene.csv \
       glycosyltransferase.txt Golgi_control_layer.txt
Rscript run_dndscv.R
python build_source_data.py exome_gene.csv Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz
python assign_tumour_type.py
Rscript per_type.R
python power2.py
python fig1_exome.py ; python figS1.py
```

Versions, install commands and per-step time and memory are in
`scripts/ENVIRONMENT.md`. The R side is pinned to an exact dndscv commit and
that pin matters: the package's version string has not changed since 2019
while its code has.

### What depends on what

| Stage | Needs | Produces |
|---|---|---|
| `get_mc3.py` | network | the MAF |
| `extract_all.py` | the MAF | `mutations_all.tsv.gz` |
| `filter_hypermutators.py` | that | `exome_filt.tsv` *(included)* |
| `build_golgi_lists.py` | `GolgiSet.xlsx` | 2 class lists *(included)* |
| `build_kinase_lists.py`, `build_kinase_classes.py` | `KinHub.pdf` | 2 class lists *(included)* |
| `run_dndscv.R` | `exome_filt.tsv`, class lists | `exome_gene.csv`, `exome_class_dnds.csv`, `genemuts.csv` *(included)* |
| `assign_tumour_type.py` | `exome_filt.tsv`, `tss.tsv` | `mut_tipo.tsv` |
| `per_type.R` | `mut_tipo.tsv`, class lists | `per_type_*.csv` *(included)* |
| `build_source_data.py` | `exome_gene.csv`, **COSMIC** | `Figure1_source_data.csv`, `cgc_by_class.csv` *(included)* |
| `power2.py` | `genemuts.csv`, `exome_gene.csv` | `power_by_gene.csv` *(included)* |
| `fig1_exome.py`, `figS1.py` | the CSVs above | the two figures |

Everything marked *(included)* is in `results/`, so any stage can be entered
from the middle without running what comes before it. The only stage that
needs a file not in this deposit is `build_source_data.py`, and its outputs
are included.

**Checkpoints.**

| After | Expect |
|---|---|
| `filter_hypermutators.py` | 2,096,036 substitutions in 10,168 tumours; `exome_filt.tsv` MD5 `d5e929031280510a816232705a8c5a4d` |
| `build_kinase_lists.py` | exactly **58** RTKs — the canonical count, so any other number means the PDF parse is wrong |
| `build_kinase_classes.py` | 435 kinases, 51 RTKs, 43 genes removed |
| `run_dndscv.R` | *TP53*, *PTEN*, *FBXW7*, *KEAP1*, *NRAS* at the top of `exome_gene.csv` by `pallsubs_cv`; 145 genes at q < 0.05 among the 19,155 with at least one substitution |
| `build_source_data.py` | **0/230, 0/269, 13/434, 6/51** |

## Known pitfalls

**cBioPortal is not a substitute for the MC3 MAF.** Its API returns no
synonymous mutations. Those are the denominator of dN/dS, so the analysis
silently yields values of 30 to 50 instead of ~1.

**MC3 writes chromosome X as `23`.** dNdScv expects `X`. Uncorrected, 5,937
mutations vanish without warning, including those in *ATP7A*, *SLC35A2*,
*C1GALT1C1* and *PIGA*.

**The background model must be fitted exome-wide.** dNdScv learns the neutral
substitution rate from the genes it is given. Restricted to the 964 genes of
interest, it gave apparent significant *negative* selection on three genes
(*CHST4*, *CHST6*, *GOPC*) and changed the class-level estimate of the Golgi
control layer from 1.023 (n.s.) to 1.043 (p = 0.014). Use `gene_list = NULL`,
then `genesetdnds`.

## Summary of results

Of the 145 genes under positive selection genome-wide, 6 are receptor tyrosine
kinases (RTK) and 13 are protein kinases. None of the 230 glycosylation enzymes
and none of the 269 Golgi trafficking/control genes is among them.

| Class | dN/dS (95% CI) | p |
|---|---|---|
| Receptor tyrosine kinases | 1.135 (1.092–1.180) | 1.0 × 10⁻¹⁰ |
| Protein kinases | 1.075 (1.056–1.095) | 1.8 × 10⁻¹⁵ |
| Golgi control layer | 1.023 (0.992–1.055) | 0.154 |
| Glycosylation enzymes | 1.031 (0.999–1.063) | 0.055 |

The Golgi control layer and the glycosylation enzymes are not significant at
p < 0.05. For a gene of typical size and mutational burden, a 2.2-fold missense
excess would have been detected; the highest observed among the 230 is 1.95
(*B3GNT1*, q = 0.43). Selection of 1.2–1.5-fold is not excluded. Analysed
within single tumour types, protein kinases or RTKs reach class-level
significance in 9 of 11 types and gene-level significance in 10 of 11; neither
Golgi class reaches significance in any type.

In the COSMIC Census, 24/51 RTKs and 59/434 kinases are listed, against 3/230
glycosylation enzymes and 7/269 Golgi genes. The recorded mechanisms for those
ten are translocation partner (*CLTC*, *CLTCL1*, *GOLGA5*, *GOPC*, *TRIP11*,
*VTI1A*, *CHST11*), recessive tumour suppressor (*EXT1*, *EXT2*) and
amplification (*GOLPH3*).

dN/dS does not capture gene dosage: a gene that is amplified but not
point-mutated returns dN/dS ≈ 1 (*GOLPH3*: 1.18, q = 0.885).

## Citation

If you use the code or the derived data, please cite both:

- **Review:** M. L. Maresca, R. Rizzo, F. Russo, S. Parashuraman and D. Russo **Tuning, not switching: the Golgi control layer as an analogue determinant of chemoresistance.**
*Bioscience Reports*, invited review (submitted 2026). https://doi.org/10.XXXX/XXXXX
- **This repository:** M. L. Maresca, R. Rizzo, F. Russo, S. Parashuraman and D. Russo (2026). *golgiproteins-vs-kinases-analysis-tcga*.
  Zenodo. https://doi.org/10.5281/zenodo.NNNNNNN
