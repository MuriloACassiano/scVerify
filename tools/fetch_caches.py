#!/usr/bin/env python
"""Build the scverify authority cache directory.

Sources (offline-first: everything is cached locally after one fetch):
  GENCODE  gencode.v44.annotation.gz  -> symbols, ENSG set, ENSG->symbol,
             ENSG->biotype, symbol->biotypes JSON files
  NCBI     Homo_sapiens.gene_info.gz  -> ncbi_gene_mapping.json
             (current symbol + every alias -> current symbol)
  HGNC     rest.genenames.org         -> hgnc_mapping.json
             (name -> current symbol or None; opt-in, network, per-symbol)

Usage:
  python tools/fetch_caches.py --only gencode
  python tools/fetch_caches.py --only ncbi
  python tools/fetch_caches.py --only hgnc --hgnc-symbols-file symbols.txt
  python tools/fetch_caches.py            # gencode + ncbi (hgnc is opt-in)
  python tools/fetch_caches.py --gtf /path/to/local/gencode.v44.annotation.gz

Downloads are kept in the cache dir as gencode.v44.annotation.gtf.gz and
Homo_sapiens.gene_info.gz so parsers can be re-run offline.
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.request

GENCODE_GTF_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/GCODE_human/"
                   "release_44/gencode.v44.annotation.gz")
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/genbank/gene/Homo_sapiens.gene_info.gz"
HGNC_API = "https://rest.genenames.org"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DEST = os.path.join(ROOT, "cache")

GTF_NAME = "gencode.v44.annotation.gtf.gz"
NCBI_NAME = "Homo_sapiens.gene_info.gz"


def http_get(url, dest, timeout=120):
    print(f"Fetching {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "scverify-cache-fetch"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as fh:
        total = r.headers.get("Content-Length")
        total = int(total) if total else None
        done = 0
        chunk = 1 << 20
        while True:
            b = r.read(chunk)
            if not b:
                break
            fh.write(b)
            done += len(b)
            if total:
                print(f"  {done / 1048576:.1f} / {total / 1048576:.1f} MB", end="\r")
    print()
    print(f"Saved {dest}")


def open_maybe_gz(path, mode="rt"):
    import gzip
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_gencode_gtf(src_path, dest_dir):
    """Parse a GENCODE annotation GTF (plain or .gz) into the five JSON caches."""
    symbols, ensgs = set(), set()
    ensg_to_symbol, ensg_to_biotype, sym_biotypes = {}, {}, {}
    n_gene = 0
    with open_maybe_gz(src_path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(ATTR_RE.findall(f[8]))
            gid = attrs.get("gene_id", "")
            name = attrs.get("gene_name", "")
            gtype = attrs.get("gene_type", "")
            if not gid.startswith("ENSG"):
                continue
            n_gene += 1
            base = gid.split(".")[0]  # canonical caches are version-less
            ensgs.add(base)
            ensg_to_biotype[base] = gtype
            if name and name != "-":
                symbols.add(name)
                ensg_to_symbol[base] = name
                sym_biotypes.setdefault(name, set()).add(gtype)
    out = {
        "gencode_v44_symbols.json": sorted(symbols),
        "gencode_v44_ensg.json": sorted(ensgs),
        "gencode_v44_ensg_to_symbol.json": ensg_to_symbol,
        "gencode_v44_ensg_to_biotype.json": ensg_to_biotype,
        "gencode_v44_symbol_to_biotypes.json":
            {k: sorted(v) for k, v in sorted(sym_biotypes.items())},
    }
    for fn, obj in out.items():
        p = os.path.join(dest_dir, fn)
        with open(p, "w") as fh:
            json.dump(obj, fh)
        print(f"Wrote {p}")
    print(f"GENCODE parse: {n_gene} gene rows, {len(symbols)} symbols, "
          f"{len(ensgs)} ENSG IDs")
    return out


def parse_ncbi_gene_info(src_path, dest_dir):
    """Parse NCBI Homo_sapiens gene_info TSV (.gz or plain) into the alias map.

    Columns (header row, '#' prefixed): col 2 = Symbol (current), col 4 =
    Synonyms (pipe-delimited aliases). Maps the symbol and every synonym to the
    current symbol.
    """
    mapping = {}
    with open_maybe_gz(src_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            symbol = parts[2].strip()
            if not symbol or symbol == "-":
                continue
            mapping[symbol] = symbol
            for a in parts[4].split("|"):
                a = a.strip()
                if a and a != "-":
                    mapping[a] = symbol
    p = os.path.join(dest_dir, "ncbi_gene_mapping.json")
    with open(p, "w") as fh:
        json.dump(mapping, fh)
    print(f"Wrote {p} ({len(mapping)} entries)")
    return mapping


def refresh_hgnc(symbols, dest_dir, sleep=0.15, timeout=15):
    """Query HGNC REST for each symbol; cache name -> current symbol or None."""
    import requests
    cache_path = os.path.join(dest_dir, "hgnc_mapping.json")
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            cache = json.load(fh)
    missing = [s for s in symbols if s not in cache]
    print(f"HGNC refresh: {len(missing)} symbols missing from cache "
          f"({len(cache)} cached)")
    for i, sym in enumerate(missing):
        try:
            resp = requests.get(f"{HGNC_API}/search/{sym}",
                                headers={"Accept": "application/json"},
                                timeout=timeout)
            if resp.status_code == 200:
                data = resp.json()
                if data["response"]["numFound"] > 0:
                    cache[sym] = data["response"]["docs"][0].get("symbol", sym)
                else:
                    cache[sym] = None
            else:
                cache[sym] = None
        except Exception as e:
            print(f"  {sym}: request failed ({e}); cached as unknown-for-now")
            cache[sym] = None
        if (i + 1) % 50 == 0:
            with open(cache_path, "w") as f:
                json.dump(cache, f, indent=2)
            print(f"  ... {i + 1}/{len(missing)}")
        time.sleep(sleep)
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2)
    print(f"HGNC cache saved: {cache_path} ({len(cache)} entries)")
    return cache


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", default=DEFAULT_DEST, help="Cache directory.")
    ap.add_argument("--gtf", default=GENCODE_GTF_URL,
                   help="GENCODE v44 GTF: URL or local path (plain or .gz).")
    ap.add_argument("--gene-info", default=NCBI_GENE_INFO_URL,
                   help="NCBI gene_info: URL or local path (.gz or plain TSV).")
    ap.add_argument("--hgnc-symbols-file", default=None,
                   help="One symbol per line; queries HGNC REST for each "
                        "uncached symbol (opt-in network step).")
    ap.add_argument("--only", default=None,
                   help="Restrict to: gencode, ncbi, or hgnc.")
    args = ap.parse_args()

    os.makedirs(args.dest, exist_ok=True)
    which = set(filter(None, [args.only])) if args.only else {"gencode", "ncbi"}

    if "gencode" in which:
        src = args.gtf
        local = os.path.join(args.dest, GTF_NAME)
        if src.startswith("http"):
            if not os.path.exists(local):
                http_get(src, local)
            src = local
        parse_gencode_gtf(src, args.dest)
    if "ncbi" in which:
        src = args.gene_info
        local = os.path.join(args.dest, NCBI_NAME)
        if src.startswith("http"):
            if not os.path.exists(local):
                http_get(src, local)
            src = local
        parse_ncbi_gene_info(src, args.dest)
    if "hgnc" in which:
        if not args.hgnc_symbols_file:
            print("--only hgnc requires --hgnc-symbols-file", file=sys.stderr)
            sys.exit(1)
        with open(args.hgnc_symbols_file) as fh:
            syms = [ln.strip() for ln in fh if ln.strip()]
        refresh_hgnc(syms, args.dest)


if __name__ == "__main__":
    main()
