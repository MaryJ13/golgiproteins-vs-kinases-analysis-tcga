# Licence

Two licences, because this deposit holds two different kinds of thing.

## Code — MIT

Everything in `scripts/` (Python, R and the helper scripts) and
`make_workdir.py`.

> Copyright (c) 2026 Domenico Russo, Consiglio Nazionale delle Ricerche
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in
> all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

## Data and documentation — CC BY 4.0

`README.md`, `METHODS.md`, everything in `results/` and `supplementary/`,
and `GolgiSet.xlsx`.

Creative Commons Attribution 4.0 International:
<https://creativecommons.org/licenses/by/4.0/>. Use it for anything, including
commercially, provided you credit the source.

---

## Third-party material, and why it is compatible

Nothing here relicenses anyone else's work. This section records what came
from elsewhere and on what terms.

**`source_data/exome_filt.tsv.gz`** — derived from the TCGA MC3 public call set,
which is open-access TCGA data. The NIH Genomic Data Sharing policy places no
restriction on redistributing open-access TCGA data; the expectation is
citation, which the README and the paper provide. Nothing in it identifies a
participant: it carries anonymised tumour barcodes, genomic coordinates and
alleles.

**`source_data/KinHub.pdf`** — the Manning kinome table as published at
kinhub.org, saved on 9 November 2026 and included unmodified so the parse can
be checked against the exact table used. The underlying content is the gene
classification of Manning *et al.*, *Science* 2002; the facts it states —
which gene belongs to which kinase group and family — are not themselves
subject to copyright, and `source_data/kinhub_parsed.csv` contains only those
facts. If you would rather not redistribute the PDF, delete it: everything
downstream runs from `kinhub_parsed.csv`, and `build_kinase_lists.py` is the
only script that needs the PDF.

**`source_data/tss.tsv`** — the TCGA tissue source site codes, a factual
lookup table published by the GDC.

**dndscv** — GPL-3, by Inigo Martincorena at the Wellcome Sanger Institute.
It is **not** redistributed here; `ENVIRONMENT.md` says how to install it from
its own repository at a pinned commit. GPL-3 therefore does not propagate to
this deposit: calling a GPL program from a separately distributed script does
not make that script a derivative work, and no dndscv code is included.

**COSMIC Cancer Gene Census** — not redistributed, by licence. Only derived
annotation for the 984 genes in the gene classes appears, in
`results/Figure1_source_data.csv` and `results/cgc_by_class.csv`: the tier,
role and mutation-type codes for those genes, the same values published in the
paper's supplementary tables. The database itself must be obtained from COSMIC
under its academic licence.

---

## What to select on Zenodo

Zenodo takes a single licence per record. Choose **CC BY 4.0**, which matches
the bulk of the deposit; the MIT terms above continue to apply to `scripts/`
and are the more permissive of the two, so nothing is lost.

If your institution requires a different arrangement, this file is the place
to change it, and the choice in the Zenodo form should be changed to match.
