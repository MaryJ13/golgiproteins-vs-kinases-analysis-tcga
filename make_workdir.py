# -*- coding: utf-8 -*-
"""
make_workdir.py — builds the flat working folder the scripts expect.

    python make_workdir.py            # creates ./work
    python make_workdir.py somewhere  # creates ./somewhere

The scripts in `scripts/` read and write the current directory. That keeps
each one short and free of path handling, but it means they cannot be run
from this folder layout directly. This script assembles the flat folder for
you, and decompresses `exome_filt.tsv.gz` so the analysis is ready to run:

    python make_workdir.py
    cd work
    python check_inputs.py

It copies rather than moves, so this deposit is left untouched and you can
delete the working folder and start again at any time. Run it twice and it
refreshes the scripts without touching anything you have produced.

Standard library only. Works on Linux, macOS and Windows.
"""
import gzip
import os
import shutil
import sys

QUI = os.path.dirname(os.path.abspath(__file__))
CARTELLE = ["scripts", "source_data", "results"]


def main():
    dest = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "work")
    if os.path.abspath(QUI) == dest:
        sys.exit("Choose a directory other than the deposit itself.")
    os.makedirs(dest, exist_ok=True)

    copiati = saltati = 0
    for c in CARTELLE:
        src = os.path.join(QUI, c)
        if not os.path.isdir(src):
            sys.exit("Cannot find %s. Run this script from inside the deposit." % src)
        for nome in sorted(os.listdir(src)):
            a = os.path.join(src, nome)
            if not os.path.isfile(a):
                continue
            b = os.path.join(dest, nome)
            if os.path.exists(b) and os.path.getmtime(b) >= os.path.getmtime(a):
                saltati += 1
                continue
            shutil.copy2(a, b)
            copiati += 1

    print("%d files copied into %s%s"
          % (copiati, dest, (", %d already current" % saltati) if saltati else ""))

    # decompress what the pipeline reads uncompressed
    for nome in ("exome_filt.tsv.gz", "per_type_genes.csv.gz"):
        gz = os.path.join(dest, nome)
        aperto = gz[:-3]
        if os.path.exists(gz) and not os.path.exists(aperto):
            with gzip.open(gz, "rb") as fi, open(aperto, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
            print("decompressed %s (%.0f MB)"
                  % (os.path.basename(aperto), os.path.getsize(aperto) / 1e6))

    print("\nnext:")
    print("  cd %s" % os.path.relpath(dest))
    print("  python check_inputs.py")


main()
