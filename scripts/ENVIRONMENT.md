# Environment

The exact versions the shipped results were produced with, and how to install
them.

## Python 3.13.18

| Package | Version | Used by |
|---|---|---|
| pandas | 3.0.5 | all the Python scripts |
| numpy | 2.5.3 | `power2.py`, `figS1.py` |
| scipy | 1.18.1 | `fig1_exome.py` (Jeffreys intervals) |
| matplotlib | 3.11.2 | `fig1_exome.py`, `figS1.py` |
| statsmodels | 0.15.0 | `power2.py`, `figS1.py` (Benjamini–Hochberg) |
| openpyxl | 3.1.5 | `build_golgi_lists.py`, to read `GolgiSet.xlsx` |
| pdfplumber | 0.11.10 | `build_kinase_lists.py`, to parse `KinHub.pdf` |
| requests | 2.34.2 | `get_mc3.py` only |

```bash
pip install "pandas==3.0.5" "numpy==2.5.3" "scipy==1.18.1" \
            "matplotlib==3.11.2" "statsmodels==0.15.0" \
            "openpyxl==3.1.5" "pdfplumber==0.11.10" "requests==2.34.2"
```

None of this is version-sensitive. Any pandas 2.x or 3.x, matplotlib 3.8+ and
statsmodels 0.13+ give the same numbers and a visually identical figure; the
Python stages were re-run on a separately built environment and reproduced
`Figure1_source_data.csv` and `power_by_gene.csv` exactly. `check_inputs.py`
and `assign_tumour_type.py` use only the standard library.

## R 4.3.3

| Package | Version | Source |
|---|---|---|
| dndscv | 0.0.1.0, commit `43c5e2f` | GitHub, `im3sanger/dndscv` |
| MASS | 7.3-60.0.1 | CRAN |
| seqinr | 4.2-36 | CRAN |
| Biostrings | 2.70.2 | Bioconductor 3.18 |
| GenomicRanges | 1.54.1 | Bioconductor 3.18 |

```r
install.packages(c("remotes", "BiocManager", "MASS", "seqinr"))
BiocManager::install(version = "3.18")
BiocManager::install(c("Biostrings", "GenomicRanges"))

# pinned to the exact commit the results were produced with
remotes::install_github("im3sanger/dndscv",
                        ref = "43c5e2f1beaaf3ac26df9c046b281ceb0edbfff7")
```

**This part is version-sensitive, and the pin matters.** dndscv's version
string has been `0.0.1.0` since 2019 while the code has continued to change —
the commit above is dated 25 August 2026 and adds a `tryCatch` around the
likelihood-ratio p value in `globaldnds`, which is one of the quantities
Figure 1D reports. Installing the package without `ref` gives whatever is
current and may not reproduce these numbers.

The package also ships its own coding-sequence reference,
`RefCDS_human_GRCh37_dNdScv.0.1.0.rda`, built on GRCh37 with the gene symbols
current at the time. The gene-list scripts map modern symbols back onto that
reference by hand, through the `ALIAS` and `TO_GRCH37` tables in
`build_golgi_lists.py` and `build_kinase_classes.py`. A different RefCDS build
would need those tables revisited and would shift the gene counts.

## R packages for the figure scripts (optional)

Only `fig1_exome.R` and `figS1.R` need these; the Python figure scripts do not.

| Package | Source |
|---|---|
| ggplot2 | CRAN |
| patchwork | CRAN |

```r
install.packages(c("ggplot2", "patchwork"))
```

## Resources

| Step | Time | Peak memory | Disk |
|---|---|---|---|
| `get_mc3.py` | 10–40 min | — | 700 MB |
| `extract_all.py` | ~15 min | < 1 GB | 40 MB |
| `filter_hypermutators.py` | ~2 min | ~1 GB | 66 MB |
| `assign_tumour_type.py` | ~1 min | < 100 MB | 125 MB |
| `run_dndscv.R` | ~50 min | ~6 GB | 4 GB for `dndsout.rds` |
| `per_type.R` | ~90 min | ~4 GB | small |
| `build_source_data.py`, `power2.py` | seconds | small | small |
| `fig1_exome.py`, `figS1.py` | seconds | small | small |

8 GB of RAM is enough; 6 GB is marginal.

`run_dndscv.R` writes `dndsout.rds`, about 4 GB, so that a mistake in the
class step does not mean paying 50 minutes for the main run again. Delete it
when you are done.
