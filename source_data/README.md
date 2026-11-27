# Source data

Six inputs in total. **Four are in this folder. Two are not, and are linked
instead** — one because it is 700 MB and already permanent at its source, the
other because its licence does not allow redistribution.

Run `python check_inputs.py` in your working folder at any point: it verifies
what is present against the checksums below and prints which parts of the
pipeline are available to you.

---

## In this folder

### `GolgiSet.xlsx` — the curated Golgi gene set
MD5 `454daef9da6884b509738f6723bfefc1` · 80 KB

507 rows, curated in-house (Russo lab, CNR Naples). **This file exists nowhere
else**; without it the gene classes cannot be rebuilt. The column `Function`
takes one of four values and the split into the two Golgi classes is made on
it:

| `Function` | rows | class |
|---|---|---|
| Glycosylation enzyme | 234 | glycosylation enzymes |
| Membrane trafficking component | 198 | Golgi control layer |
| Channels and transporters | 40 | Golgi control layer |
| lipid trafficking and metabolism | 35 | Golgi control layer |

### `KinHub.pdf` — the human kinome
MD5 `3dd2b0fe28d8b6f98d0ba4c3bb4ceb28` · 2.2 MB

The Manning kinome as published at <http://www.kinhub.org/kinases.html>, saved
as PDF on 9 November 2026. Kept as a file rather than as a link so the parse
can be checked against *this* table rather than a later revision of the site.
519 kinases with group and family. `build_kinase_lists.py` parses it, and the
parse is correct if and only if it reports exactly 58 receptor tyrosine
kinases.

### `kinhub_parsed.csv`
That PDF already parsed, so the parsing step can be skipped or checked.

### `tss.tsv`
TCGA tissue source site codes mapped to tumour types, from the GDC code table
at <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables>. Used only
by `assign_tumour_type.py`. All 10,168 tumours map with no unrecognised codes.

### `exome_filt.tsv.gz` — the mutations, already filtered
MD5 of the uncompressed file `d5e929031280510a816232705a8c5a4d` · 14 MB
gzipped, 66 MB open

It is the MC3 call set after extraction and filtering: 2,096,036 single-base
substitutions in 10,168 tumours, in the five columns dNdScv consumes
(`sampleID`, `chr`, `pos`, `ref`, `mut`). It is derived from open-access TCGA
data (see `LICENSE.md`); at 14 MB it saves a 700 MB download and about 35 
minutes of processing.

With this file you can go straight to `Rscript run_dndscv.R` and reproduce
everything downstream. Download the MAF only if you want to verify the
extraction itself.

```bash
gunzip exome_filt.tsv.gz
```

---

## Linked, not included

### 1. `mc3.v0.2.8.PUBLIC.maf.gz` — the full call set
About 700 MB

```
https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
```

`python get_mc3.py` fetches exactly that URL and resumes if interrupted. The
dataset's landing page is
<https://gdc.cancer.gov/about-data/publications/mc3-2017>.

Not included because it is large, already public, and permanent at the GDC,
where it is versioned and citable. Shipping a copy would add a second version
that can drift from the canonical one.

**Needed only to rebuild `exome_filt.tsv` from scratch.** That file is already
here, so for every other purpose this download is optional.

**Do not substitute cBioPortal.** its API returns no synonymous mutations, which are
the denominator of dN/dS, and an analysis built on it gives dN/dS values of 30
to 50 instead of about 1.

After extraction, expect 3,284,873 substitutions in 10,294 tumours, which
`filter_hypermutators.py` reduces to 2,096,036 in 10,168. A different count 
indicates a problem in extraction or filtering.

### 2. COSMIC Cancer Gene Census, release v104 (19 May 2026)

```
https://cancer.sanger.ac.uk/cosmic/download/cosmic
```

Free for academic use behind a login, and **not redistributable** — this is a
licence constraint, not a size one, so it is not included. Take release v104 and 
the file `Cosmic_CancerGeneCensus_Tsv_v104_GRCh37.tar`; untar it and leave the
`.tsv.gz` in the working folder.

| File | MD5 of the `.tsv.gz` |
|---|---|
| `Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz` | `5b07d982cf63851c56bea1728557c9c1` |
| `Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz` | `84b36425c805f1b49c1083991c3eb266` |

Both releases were checked and give identical answers for every gene reported
here, because the analysis reads only `GENE_SYMBOL`, `TIER`, `ROLE_IN_CANCER`
and `MUTATION_TYPES`, never a coordinate. **The GRCh37 file alone is
sufficient**; the GRCh38 one was used only as a cross-check.

**Needed only by `build_source_data.py`.** Its outputs —
`Figure1_source_data.csv` and `cgc_by_class.csv`, which carry the Census
annotation for the 984 class genes — are in `results/`, so every panel of both
figures can be produced without this download. You need it only to regenerate
that annotation from the database itself.

Release v104 is used because the release bundled with dNdScv (v81, 2017) does
not list *EXT1*, *EXT2*, *CHST11* or *GOLPH3* in the Census. In v104, *EXT1* and
*EXT2* are Tier 1, *CHST11* is Tier 2, and *GOLPH3* is listed with amplification
as its only mutation type.

---

## Reference genome

Everything is **GRCh37**. The MC3 call set is GRCh37, and the coding sequence
reference bundled with dNdScv (`RefCDS_human_GRCh37_dNdScv.0.1.0.rda`) is
GRCh37. No coordinate conversion happens anywhere in the pipeline, and none
should be introduced.
