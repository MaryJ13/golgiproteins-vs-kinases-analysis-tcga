# -*- coding: utf-8 -*-
"""
Downloads the TCGA MC3 public MAF.

MC3 is the consortium-harmonised pan-cancer call set: 33 tumour types, around
10,000 tumours, aligned to GRCh37, and it contains every variant class.
Synonymous substitutions matter here because they are the denominator of
dN/dS. The cBioPortal API returns only non-synonymous calls, which makes
dN/dS incomputable from that source and silently yields ratios of 30 to 50
instead of about 1.

One file, one URL:
    https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
    -> mc3.v0.2.8.PUBLIC.maf.gz   (about 700 MB)

    pip install requests
    python get_mc3.py

The download resumes if interrupted. The MAF is kept on disk; extract_all.py
reads it next. Delete it by hand when the pipeline has finished with it.

The landing page for the dataset is
https://gdc.cancer.gov/about-data/publications/mc3-2017
"""
import os
import sys

try:
    import requests
except ImportError:
    sys.exit("First run:  pip install requests")

URL = "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"
MAF = "mc3.v0.2.8.PUBLIC.maf.gz"


def dimensione_remota():
    try:
        r = requests.head(URL, allow_redirects=True, timeout=120)
        return int(r.headers.get("Content-Length", 0))
    except Exception:
        return 0


def scarica(dest):
    totale = dimensione_remota()
    fatto = os.path.getsize(dest) if os.path.exists(dest) else 0

    if totale and fatto == totale:
        print("%s already complete (%.0f MB)" % (dest, totale / 1e6))
        return
    if fatto:
        print("resuming from %.0f MB" % (fatto / 1e6))

    testa = {"Range": "bytes=%d-" % fatto} if fatto else {}
    with requests.get(URL, headers=testa, stream=True, timeout=600) as r:
        if fatto and r.status_code == 200:
            # the server ignored the range request: start again
            fatto = 0
        elif fatto and r.status_code != 206:
            r.raise_for_status()
        else:
            r.raise_for_status()
        modo = "ab" if fatto else "wb"
        with open(dest, modo) as fh:
            for blocco in r.iter_content(chunk_size=1 << 20):
                if not blocco:
                    continue
                fh.write(blocco)
                fatto += len(blocco)
                if totale:
                    print("\r  %.0f / %.0f MB  (%.1f%%)"
                          % (fatto / 1e6, totale / 1e6, 100 * fatto / totale),
                          end="", flush=True)
                else:
                    print("\r  %.0f MB" % (fatto / 1e6), end="", flush=True)
    print()


def main():
    dest = os.path.join(os.getcwd(), MAF)
    scarica(dest)
    n = os.path.getsize(dest)
    print("\nwritten %s  (%.0f MB)" % (dest, n / 1e6))
    if n < 500e6:
        print("WARNING: the file looks too small. Expected about 700 MB.")
        print("Delete it and run this script again.")
        return
    print("\nnext: python extract_all.py")


main()
