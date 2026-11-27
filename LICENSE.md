# Licence

Two licences, because this deposit holds two different kinds of thing.

## Code — MIT

Everything in `scripts/` (Python, R and the helper scripts) and
`make_workdir.py`.

> Copyright (c) 2026 Maresca Maria Luigia, Rizzo Riccardo, Russo Francesco, Parashuraman Seetharaman and Russo Domenico, Consiglio Nazionale delle Ricerche
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
and `source_data/GolgiSet.xlsx`.

Creative Commons Attribution 4.0 International:
<https://creativecommons.org/licenses/by/4.0/>. Reuse, including commercial reuse, is permitted with attribution.

---

## Third-party material

Nothing in this deposit relicenses anyone else's work. This section records the source of each third-party item and the terms attached to it by its source.

**`source_data/exome_filt.tsv.gz`** -> derived from the TCGA MC3 public call set (mc3.v0.2.8.PUBLIC.maf.gz, open-access data from the NCI Genomic Data Commons; Ellrott et al., Cell Systems 2018). It contains anonymised tumour barcodes, genomic coordinates and alleles.

**`source_data/KinHub.pdf`** -> the kinome table published at kinhub.org (http://www.kinhub.org/kinases.html), saved on 9 November 2026 and included unmodified so that the parse can be checked against the table used. The classification is that of Manning et al., Science 2002. source_data/kinhub_parsed.csv contains the gene symbols with their group and family. build_kinase_lists.py is the only script that reads the PDF; everything downstream runs from kinhub_parsed.csv.

**`source_data/tss.tsv`** -> the TCGA tissue source site codes, from the GDC code tables (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables).

**dndscv** -> GPL-3 licence (Martincorena et al., Cell 2017). Not included in this deposit; scripts/ENVIRONMENT.md gives the installation from its own repository at the commit used.

**COSMIC Cancer Gene Census** -> not included; it is distributed by COSMIC under its own licence (https://cancer.sanger.ac.uk/cosmic). results/Figure1_source_data.csv and results/cgc_by_class.csv contain derived annotation (tier, role in cancer, mutation-type codes) for the 984 genes in the four gene classes.

---
