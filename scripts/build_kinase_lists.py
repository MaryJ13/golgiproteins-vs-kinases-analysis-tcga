#!/usr/bin/env python3
"""
Build genes_kinases.txt and genes_rtk.txt from the KinHub table.

KinHub (http://www.kinhub.org/kinases.html) is the Manning kinome with HGNC
symbols attached. The site is HTTP-only and awkward to scrape, so this script
reads the saved PDF of that page instead, which makes the lists reproducible
offline.

    pip install pdfplumber        # or just have pdftotext on the PATH
    python build_kinase_lists.py KinHub.pdf

KinHub's symbols are mostly current but the table predates several HGNC
renames, and two rows have an empty HGNC cell. ALIASES below fixes both. Any
symbol cBioPortal cannot map is reported by fetch_cbioportal.py, so check that
report rather than assuming this list is clean.

Receptor tyrosine kinases are the TK group minus the nine non-receptor
families. That rule reproduces the canonical set of 58 human RTKs exactly,
which is the check that the parse worked.

Receptor serine/threonine kinases (TGFBR, BMPR, ACVR) and the guanylate
cyclase receptors sit in TKL and RGC, not TK, so they stay in the kinase
class. That is deliberate: the figure contrasts tyrosine-kinase receptors with
the Golgi classes.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

GROUPS = {"TK", "TKL", "AGC", "CAMK", "CMGC", "CK1", "STE", "Atypical",
          "RGC", "Other"}

# TK families whose members are cytoplasmic, not receptors.
NON_RECEPTOR_TK = {"Abl", "Ack", "Csk", "FAK", "Fer", "Jak", "Src", "Syk", "Tec"}

# KinHub symbol -> current HGNC symbol.
ALIASES = {
    "ADRBK1": "GRK2",      # beta-adrenergic receptor kinase 1
    "ADRBK2": "GRK3",
    "COL4A3BP": "CERT1",   # ceramide transfer protein; also in the GolgiSet
    "GSG2": "HASPIN",
    "ADCK3": "COQ8A",
    "ADCK4": "COQ8B",
    "ICK": "CILK1",
    "PAK7": "PAK5",
    "ZAK": "MAP3K20",
}

# Rows whose HGNC cell is empty in the PDF, keyed by the xName column.
BY_XNAME = {
    "ZAK": ("MAP3K20", "TKL", "MLK"),
    "SgK494": ("RSKR", "AGC", "RSKR"),
}

JUNK = {"nan", "Obscurin", "Name", "Gene", "Kinase", "Nomenclature"}


def pdf_lines(path):
    """Layout-preserving text, via pdfplumber if present, else pdftotext."""
    try:
        import pdfplumber
        out = []
        with pdfplumber.open(path) as pdf:
            for page in pdf.pages:
                out += (page.extract_text(layout=True) or "").split("\n")
        if any(l.strip() for l in out):
            return out
    except ImportError:
        pass
    try:
        txt = subprocess.run(["pdftotext", "-layout", str(path), "-"],
                             capture_output=True, text=True, check=True).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        sys.exit(f"Could not read the PDF. Install pdfplumber (pip install "
                 f"pdfplumber) or poppler-utils for pdftotext. ({e})")
    return txt.split("\n")


def parse(lines):
    """Token-based, so it does not depend on how the extractor spaces columns.

    Every data row begins with three single-word identifiers, xName, the
    Manning name and the HGNC symbol, followed by the free-text kinase name
    and then the Group token. Splitting on whitespace and taking the first
    three tokens therefore survives both pdftotext and pdfplumber, which
    disagree about how many spaces sit between columns.
    """
    rows = {}
    for line in lines:
        if not line.strip() or "KinHub" in line or "Manning et al" in line:
            continue
        t = line.split()
        if len(t) < 4:
            continue
        gi = next((i for i, p in enumerate(t) if p in GROUPS and i >= 3), None)
        if gi is None:
            continue

        xname, sym = t[0], t[2]
        group = t[gi]
        family = t[gi + 1] if len(t) > gi + 1 else ""

        if xname in BY_XNAME:          # HGNC cell is blank in the PDF
            sym, group, family = BY_XNAME[xname]
        if sym in JUNK or not re.fullmatch(r"[A-Z][A-Z0-9\-]{1,9}", sym):
            continue
        sym = ALIASES.get(sym, sym)
        rows.setdefault(sym, {"symbol": sym, "xname": xname,
                              "group": group, "family": family})
    return list(rows.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", nargs="?", default="KinHub.pdf")
    ap.add_argument("--genedir", default=".")
    args = ap.parse_args()

    if not Path(args.pdf).exists():
        sys.exit(f"{args.pdf} not found. Save http://www.kinhub.org/kinases.html "
                 f"as a PDF and pass its path.")

    rows = parse(pdf_lines(args.pdf))
    if len(rows) < 400:
        sys.exit(f"Only {len(rows)} kinases parsed; the PDF layout is not what "
                 f"this script expects. Inspect kinhub_parsed.csv.")

    tk = [r for r in rows if r["group"] == "TK"]
    rtk = sorted(r["symbol"] for r in tk
                 if r["family"] not in NON_RECEPTOR_TK)
    kinases = sorted(r["symbol"] for r in rows if r["symbol"] not in set(rtk))

    # Golgi genes win: a GolgiSet gene must not also serve as a control.
    golgi = set()
    for f in ("genes_glycosyltransferases.txt", "genes_golgi_control_layer.txt"):
        p = Path(args.genedir) / f
        if p.exists():
            golgi |= {l.strip() for l in p.read_text().splitlines()
                      if l.strip() and not l.startswith("#")}
    clash = sorted((set(kinases) | set(rtk)) & golgi)
    if clash:
        print(f"{len(clash)} gene(s) are in both the GolgiSet and the kinome: "
              f"{', '.join(clash)}")
        print("Dropped from the kinase classes so the four classes stay disjoint.")
        kinases = [g for g in kinases if g not in golgi]
        rtk = [g for g in rtk if g not in golgi]

    hdr_k = (f"# Protein kinases from KinHub ({args.pdf}), receptor tyrosine "
             f"kinases excluded.\n"
             f"# Receptor Ser/Thr kinases (TGFBR, BMPR, ACVR) and guanylate "
             f"cyclase receptors are here,\n"
             f"# not in the RTK list: they are TKL and RGC, not TK.\n")
    hdr_r = (f"# Receptor tyrosine kinases: KinHub TK group minus the "
             f"non-receptor families\n"
             f"# ({', '.join(sorted(NON_RECEPTOR_TK))}).\n")
    Path(args.genedir, "genes_kinases.txt").write_text(
        hdr_k + "\n".join(kinases) + "\n")
    Path(args.genedir, "genes_rtk.txt").write_text(
        hdr_r + "\n".join(rtk) + "\n")

    import csv
    with open(Path(args.genedir, "kinhub_parsed.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, ["symbol", "xname", "group", "family"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["group"], r["family"],
                                                r["symbol"])))

    print(f"\n{len(rows)} unique kinases parsed")
    print(f"  genes_rtk.txt      {len(rtk)}")
    print(f"  genes_kinases.txt  {len(kinases)}")
    print(f"  kinhub_parsed.csv  full table with group and family")
    if len(rtk) != 58:
        print(f"\n! {len(rtk)} RTKs, expected 58. Check the Family column in "
              f"kinhub_parsed.csv against NON_RECEPTOR_TK.")
    else:
        print("\nRTK count is 58, the canonical human set: the parse is sound.")


if __name__ == "__main__":
    main()
