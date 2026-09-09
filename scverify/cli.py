"""Command-line interface for scverify.

READ-ONLY. Never modifies, filters, or deletes data. Exit code is 1 on FAIL,
2 on WARN with --strict, else 0.

Usage:
  python -m scverify GSE262381
  python -m scverify --all
  python -m scverify --all --out tmp/v
  python -m scverify path/to/dataset.h5ad --no-distro --brief
  python -m scverify --all --strict
  python -m scverify GSE178341 --cache-dir /path/to/cache
"""

import argparse
import os
import sys

from scverify.verification import (
    CACHE_DIR,
    discover_feature_reference,
    load_adt_panel,
    load_gencode_ensg,
    load_gencode_maps,
    load_gencode_symbols,
    load_hgnc_cache,
    load_ncbi_cache,
    refresh_hgnc_cache,
    test_title,
    verify_dataset,
)


def main():
    parser = argparse.ArgumentParser(
        description="Read-only ADT-in-GEX verification: reference tests + "
                    "empirical ADT-distribution tests.")
    parser.add_argument("dataset", nargs="?",
                       help="Dataset directory label (with --base-dir layout) or a "
                            "direct path to a .h5ad file.")
    parser.add_argument("--all", action="store_true",
                        help="Verify every GSE* directory under --base-dir.")
    parser.add_argument("--base-dir", default="file_database",
                        help="Root containing one directory per dataset "
                             "(default: ./file_database).")
    parser.add_argument("--cache-dir", default=CACHE_DIR,
                        help="Authority cache directory (default: <project>/cache).")
    parser.add_argument("--out", default=None,
                        help="Write per-dataset standard + test tables to CSV "
                             "(<out>_<DS>_table.csv, <out>_<DS>_tests.csv).")
    parser.add_argument("--no-distro", action="store_true",
                        help="Skip full per-feature distribution statistics "
                             "(DIST tests reduce to ADT-name counts).")
    parser.add_argument("--chunk-size", type=int, default=500,
                        help="Rows per chunk for distribution stats (default 500).")
    parser.add_argument("--verbose", action="store_true",
                        help="Print full tables (default is a compact digest; "
                             "full tables always go to --out CSVs).")
    parser.add_argument("--brief", action="store_true",
                        help="Ultra-brief output: only OVERALL, placements, and "
                             "WARN/FAIL one-liners.")
    parser.add_argument("--strict", action="store_true",
                        help="WARN also yields a nonzero exit code.")
    parser.add_argument("--refresh-hgnc", action="store_true",
                        help="Opt-in: query HGNC REST for symbols missing from "
                             "the local cache before verifying.")
    args = parser.parse_args()

    cache_dir = os.path.abspath(args.cache_dir)
    base_dir = os.path.abspath(args.base_dir)

    gencode_symbols = load_gencode_symbols(cache_dir)
    gencode_ensg = load_gencode_ensg(cache_dir)
    gmaps = load_gencode_maps(cache_dir)
    print(f"GENCODE v44: {len(gencode_symbols)} symbols, {len(gencode_ensg)} ENSG IDs, "
          f"maps: {len(gmaps['ensg_to_symbol'])} ensg->symbol, "
          f"{len(gmaps['ensg_to_biotype'])} ensg->biotype")
    hgnc_cache = load_hgnc_cache(cache_dir)
    ncbi_cache = load_ncbi_cache(cache_dir)
    adt_panel = load_adt_panel(cache_dir)
    print(f"Caches: HGNC={len(hgnc_cache)} entries, NCBI={len(ncbi_cache)} entries "
          "(offline; use --refresh-hgnc to query HGNC)")
    if adt_panel:
        print(f"ADT panel: {adt_panel['provenance']}")
    else:
        print("ADT panel: MISSING (protein-envelope screen and sparsity check will WARN)")
    print("How to read: per dataset, OVERALL says CLEAN / PASS WITH NOTES / FAIL. "
          "Only WARN/FAIL checks print details; full tables via --verbose or --out CSVs.")

    if args.all:
        if not os.path.isdir(base_dir):
            print(f"error: --base-dir {base_dir} does not exist", file=sys.stderr)
            sys.exit(1)
        jobs = [{"label": d, "h5ad": None, "raw_dir": None}
                for d in sorted(os.listdir(base_dir))
                if d.startswith("GSE")
                and os.path.isdir(os.path.join(base_dir, d))]
    elif args.dataset:
        if args.dataset.endswith(".h5ad") and os.path.exists(args.dataset):
            # Direct file mode: verify any h5ad in place (no layout assumed).
            ap = os.path.abspath(args.dataset)
            stem = os.path.splitext(os.path.basename(ap))[0]
            sibling_raw = os.path.join(os.path.dirname(ap), "raw")
            jobs = [{"label": stem, "h5ad": ap,
                     "raw_dir": sibling_raw
                     if os.path.isdir(sibling_raw) else None}]
        else:
            jobs = [{"label": args.dataset, "h5ad": None, "raw_dir": None}]
    else:
        parser.print_help()
        sys.exit(1)

    if args.refresh_hgnc:
        import anndata as ad
        for job in jobs:
            p = job["h5ad"] or os.path.join(base_dir, job["label"],
                                            f"{job['label']}.h5ad")
            if os.path.exists(p):
                a = ad.read_h5ad(p, backed="r")
                try:
                    refresh_hgnc_cache(a.var_names.tolist(), cache_dir)
                finally:
                    a.file.close()
        hgnc_cache = load_hgnc_cache(cache_dir)

    results = {}
    for job in jobs:
        ds = job["label"]
        ref_map, ref_provenance = discover_feature_reference(
            ds, base_dir, cache_dir=cache_dir, raw_dir=job["raw_dir"])
        print(f"Feature reference for {ds}: "
              f"{len(ref_map) if ref_map else 0} rows "
              f"({ref_provenance})" if ref_map else
              f"Feature reference for {ds}: none ({ref_provenance})")
        r = verify_dataset(ds, base_dir, gencode_symbols, gencode_ensg, gmaps,
                           hgnc_cache, ncbi_cache, adt_panel,
                           chunk_size=args.chunk_size,
                           run_distro=not args.no_distro,
                           ref_map=ref_map, ref_provenance=ref_provenance,
                           h5ad_path=job["h5ad"], verbose=args.verbose,
                           brief=args.brief)
        if r:
            results[ds] = r
            if args.out:
                base, ext = os.path.splitext(args.out)
                ext = ext or ".csv"
                d = os.path.dirname(base)
                if d:
                    os.makedirs(d, exist_ok=True)
                tp = f"{base}_{ds}_table.csv"
                tt = f"{base}_{ds}_tests.csv"
                r["standard_table"].to_csv(tp, index=False)
                r["test_table"].to_csv(tt, index=False)
                print(f"  Wrote {tp} and {tt}")

    print(f"\n\n{'=' * 70}")
    print("  SUMMARY")
    print(f"{'=' * 70}")
    any_fail, any_warn = False, False
    for ds, r in results.items():
        t = r["test_table"]
        fails = t[t["result"] == "FAIL"]
        warns = t[t["result"] == "WARN"]
        any_fail |= len(fails) > 0
        any_warn |= len(warns) > 0
        status = ("FAIL:" + ",".join(test_title(t) for t in fails["test"]) if len(fails)
                  else ("WARN:" + ",".join(test_title(t) for t in warns["test"]) if len(warns) else "ALL PASS"))
        print(f"  {ds}: {r['n_cells']} cells x {r['n_genes']} genes | {status}")
        print(f"    OVERALL: {r['overall']}")
        print(f"    {r['adt_placement']}")
        print(f"    {r['mux_placement']}")

    if any_fail:
        sys.exit(1)
    if any_warn and args.strict:
        sys.exit(2)


if __name__ == "__main__":
    main()
