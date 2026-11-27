# -*- coding: utf-8 -*-
"""
check_inputs.py — what is present, and what you can therefore run.

    python check_inputs.py

Run it in the working folder before anything else. It verifies the files that
are present against their recorded checksums and prints which stages of the
pipeline are available to you. It downloads nothing and changes nothing.

There are three ways to use this deposit, and the point of this script is to
tell you which ones are open:

  1. REDRAW   the figures from the shipped results. No downloads at all.
  2. REANALYSE from the filtered mutation table, which is included.
              Skips the 700 MB download and the extraction. Needs R.
  3. REBUILD  everything from the original sources. Needs the MC3 MAF and,
              for the Cancer Gene Census annotation only, the COSMIC release.
"""
import hashlib
import os

ATTESI = {
    # file                              md5                                  what it is
    "GolgiSet.xlsx":        ("454daef9da6884b509738f6723bfefc1", "curated Golgi gene set (included)"),
    "KinHub.pdf":           ("3dd2b0fe28d8b6f98d0ba4c3bb4ceb28", "the kinome as parsed (included)"),
    "exome_filt.tsv":       ("d5e929031280510a816232705a8c5a4d", "filtered substitutions, ready for dNdScv (included, gzipped)"),
}
OPZIONALI = {
    "mc3.v0.2.8.PUBLIC.maf.gz": "the full MC3 call set (download; see source_data/README.md)",
    "Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz": "Cancer Gene Census (download under academic licence)",
}
RISULTATI = ["exome_gene.csv", "Figure1_source_data.csv", "exome_class_dnds.csv",
             "cgc_by_class.csv", "power_by_gene.csv", "per_type_class.csv",
             "genemuts.csv"]
LISTE = ["glycosyltransferase.txt", "Golgi_control_layer.txt", "kinase.txt",
         "receptor_tyrosine_kinase.txt"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for blocco in iter(lambda: fh.read(1 << 20), b""):
            h.update(blocco)
    return h.hexdigest()


def riga(stato, nome, nota=""):
    print("  [%s] %-46s %s" % (stato, nome, nota))


def main():
    print("\nINPUTS\n")
    ok = {}
    for nome, (atteso, cosa) in ATTESI.items():
        p = nome if os.path.exists(nome) else (nome + ".gz" if os.path.exists(nome + ".gz") else None)
        if p is None:
            riga(" ", nome, "MISSING -- " + cosa)
            ok[nome] = False
        elif p.endswith(".gz") and not nome.endswith(".gz"):
            riga("~", p, "present, still gzipped: run  gunzip " + p)
            ok[nome] = False
        else:
            h = md5(p)
            if h == atteso:
                riga("x", nome, "verified")
            else:
                riga("!", nome, "CHECKSUM MISMATCH -- expected %s, got %s" % (atteso, h))
            ok[nome] = (h == atteso)

    for nome, cosa in OPZIONALI.items():
        if os.path.exists(nome):
            riga("x", nome, "present")
            ok[nome] = True
        else:
            riga(" ", nome, "not present -- " + cosa)
            ok[nome] = False

    print("\nSHIPPED RESULTS\n")
    mancanti = [f for f in RISULTATI + LISTE if not os.path.exists(f)]
    if mancanti:
        riga("!", "%d of %d missing" % (len(mancanti), len(RISULTATI) + len(LISTE)),
             ", ".join(mancanti[:4]))
    else:
        riga("x", "all %d present" % (len(RISULTATI) + len(LISTE)), "")

    print("\nWHAT YOU CAN RUN\n")
    redraw = not mancanti
    reanalyse = ok.get("exome_filt.tsv", False) and not mancanti
    rebuild = ok.get("mc3.v0.2.8.PUBLIC.maf.gz", False)
    census = ok.get("Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz", False)

    print("  1. REDRAW the figures            %s" % ("YES  -> python fig1_exome.py ; python figS1.py"
                                                     if redraw else "no, shipped results incomplete"))
    print("  2. REANALYSE from exome_filt.tsv %s" % ("YES  -> Rscript run_dndscv.R  (~50 min, ~6 GB RAM)"
                                                     if reanalyse else "no, exome_filt.tsv missing or unverified"))
    print("  3. REBUILD from the MAF          %s" % ("YES  -> python extract_all.py"
                                                     if rebuild else "no, download the MC3 MAF first (python get_mc3.py)"))
    print("     rebuild Census annotation     %s" % ("YES  -> python build_source_data.py ..."
                                                     if census else "no, needs the COSMIC release; the shipped\n"
                                                     "                                      Figure1_source_data.csv and cgc_by_class.csv\n"
                                                     "                                      already carry that annotation, so every\n"
                                                     "                                      figure panel can be made without it"))
    print()


main()
