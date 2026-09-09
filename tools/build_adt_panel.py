#!/usr/bin/env python
"""Build cache/adt_reference_panel.json from a CITE-seq h5ad's obsm matrix.

The panel is the empirical reference for the DIST tests: per-antibody
moment statistics plus a p5-p95 envelope, calibrated on REAL antibody
counts, not invented thresholds.

Usage:
  python tools/build_adt_panel.py /path/to/cite_seq.h5ad
  python tools/build_adt_panel.py GSE262381.h5ad --obsm-key adt \\
      --out cache/adt_reference_panel.json
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import anndata as ad
import numpy as np
from scipy.stats import kurtosis, skew

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUT = os.path.join(ROOT, "cache", "adt_reference_panel.json")


def pct(x, q):
    return float(np.percentile(x, q))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("h5ad", help="CITE-seq h5ad with the antibody matrix in obsm.")
    ap.add_argument("--obsm-key", default="adt",
                   help="obsm key holding the antibody matrix (default: adt).")
    ap.add_argument("--out", default=DEFAULT_OUT,
                   help=f"Output JSON (default: {DEFAULT_OUT}).")
    ap.add_argument("--provenance", default=None,
                   help="Provenance string (default: auto from the file).")
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        if args.obsm_key not in adata.obsm:
            print(f"error: obsm['{args.obsm_key}'] not found; "
                  f"available: {list(adata.obsm.keys())}", file=sys.stderr)
            sys.exit(1)
        raw = adata.obsm[args.obsm_key]
        if hasattr(raw, "columns"):          # pandas DataFrame
            labels = [str(x) for x in raw.columns]
            M = np.asarray(raw.to_numpy(), dtype=float)
        elif hasattr(raw, "var_names"):      # sub-AnnData
            labels = [str(x) for x in raw.var_names]
            M = np.asarray(raw.X.toarray() if hasattr(raw.X, "toarray")
                           else raw.X, dtype=float)
        else:                                # plain array
            labels = []
            M = np.asarray(raw, dtype=float)
        n_cells, n_ab = M.shape
        if args.provenance is None:
            label = os.path.splitext(os.path.basename(args.h5ad))[0]
            args.provenance = (f"{label} obsm[{args.obsm_key}], "
                               f"{n_ab} antibodies x {n_cells} cells")
        if len(labels) != n_ab:
            labels = [f"ab{i}" for i in range(n_ab)]
    finally:
        adata.file.close()

    features = []
    for j in range(n_ab):
        col = M[:, j].astype(float)
        mean = float(col.mean())
        var = float(col.var())
        frac = float((col > 0).mean())
        with np.errstate(all="ignore"):
            sk = skew(col)
            ku = kurtosis(col, fisher=False)
            bc = (sk ** 2 + 1) / ku if ku not in (0, 0.0) and np.isfinite(ku) else np.nan
        features.append({
            "antibody": labels[j],
            "mean": mean,
            "var": var,
            "overdisp": float(var / mean) if mean > 0 else 0.0,
            "frac_expressed": frac,
            "bc": float(bc) if np.isfinite(bc) else None,
            "median": float(np.median(col)),
            "max": float(col.max()),
        })

    stats = np.array([[f["mean"], f["overdisp"], f["frac_expressed"],
                       np.nan if f["bc"] is None else f["bc"]]
                      for f in features])
    keys = ["mean", "overdisp", "frac_expressed", "bc"]
    envelope = {}
    for k, col in zip(keys, stats.T):
        col = col[np.isfinite(col)]
        if col.size:
            envelope[k] = {"p5": pct(col, 5), "p50": pct(col, 50),
                          "p95": pct(col, 95)}
        else:
            envelope[k] = {"p5": None, "p50": None, "p95": None}

    panel = {
        "provenance": args.provenance,
        "n_cells": int(n_cells),
        "n_antibodies": int(n_ab),
        "features": features,
        "envelope_p5_p95": envelope,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(panel, fh, indent=1)
    print(f"Wrote {args.out}: {n_ab} antibodies x {n_cells} cells "
          f"({args.provenance})")


if __name__ == "__main__":
    main()
