# scVerify

Read-only verification that antibody-derived tags (ADT/CITE-seq) and multiplexing
hashtags (HTO) are **not mixed into the RNA (GEX) matrix** of single-cell
`.h5ad` files — and that feature names resolve against authoritative references.

scVerify **classifies, tests, and reports only**. It never modifies, filters,
renames, or deletes data. Every finding is `PASS` / `WARN` / `FAIL` with a
stable test ID, an affected-count, and examples; every feature is retained.

## What it checks

- **REF** — feature names vs GENCODE v44 (symbols, ENSG IDs, biotypes), the
  local HGNC cache, and the local NCBI Gene cache, plus symbol↔ENSG
  cross-column consistency and feature-reference (barcode) reconciliation.
- **DIST** — per-feature count distributions vs an **empirical ADT reference
  panel** (137 real antibody columns from GSE262381 `obsm['adt']`), and
  ADT-spelled features with nonzero RNA counts (the core ADT-in-GEX check).
- **MUX** — multiplexing/hashtag tag discovery, placement, and demux
  bookkeeping consistency, including unsupervised, counts-only hashtag
  discovery.

See `docs/METHOD.md` for the full test catalogue and evidence tiers, and
`docs/MANUAL.md` / `docs/OUTPUTS.md` for the CLI and output formats.

## Installation

Python ≥ 3.10. Dependencies: numpy, pandas, scipy, anndata, h5py, requests.

Validated under both pandas 2.3 (Python 3.11, anndata 0.12) and pandas 3.0
(Python 3.14, anndata 0.13); reports are byte-identical across the two.

```bash
# conda (matches environment.yml)
conda env create -f environment.yml
conda activate scverify

# or, in any prepared environment:
pip install -e .
```

The `cache/` directory ships with the authority caches prebuilt (GENCODE v44,
NCBI Gene, HGNC, ADT panel, plus the raw GTF/gene_info sources), so the tool
runs **fully offline** out of the box. Use `tools/fetch_caches.py` to rebuild
or refresh them (see `docs/MANUAL.md`).

## Quick start

```bash
# direct file mode (any .h5ad, in place; a sibling raw/ dir is used if present)
python -m scverify /path/to/dataset.h5ad --brief

# full per-feature distribution tests (slower; use --no-distro to skip)
python -m scverify /path/to/dataset.h5ad

# base-dir layout: one subdirectory per dataset, <DS>.h5ad inside
python -m scverify --all --base-dir /path/to/file_database --no-distro --brief

# console script (after pip install -e .)
scverify /path/to/dataset.h5ad --no-distro --brief
```

Exit code: `1` if any test `FAIL`s, `2` if only `WARN`s and `--strict` is set,
otherwise `0`.

## Project layout

```
scverify/            Python package
  verification.py    core: classification, REF/DIST/MUX tests, cache loaders
  cli.py             command-line interface
  __main__.py        python -m scverify
tools/
  fetch_caches.py    rebuild GENCODE/NCBI caches (offline-capable) + HGNC refresh
  build_adt_panel.py rebuild cache/adt_reference_panel.json from a CITE-seq h5ad
cache/               authority caches (ship with the project; offline-first)
docs/                METHOD.md, MANUAL.md, OUTPUTS.md
```

## License

CC BY-NC-SA 4.0 — see `LICENSE`.
