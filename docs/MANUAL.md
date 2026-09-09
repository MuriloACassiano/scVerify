# Manual

## Invocation

```bash
python -m scverify [dataset] [options]      # module form
scverify [dataset] [options]                # console script (after pip install -e .)
```

## Arguments

| Argument | Meaning |
|----------|---------|
| `dataset` | Optional positional: a single `.h5ad` file. Omit when using `--all`. |
| `--all` | Discover every dataset under `--base-dir` and verify each. |
| `--base-dir DIR` | Base directory for dataset discovery (default: `file_database`). Each subdirectory containing a `.h5ad` is one dataset; a sibling `raw/` directory, if present, supplies the feature reference. |
| `--cache-dir DIR` | Authority cache directory (default: the project `cache/`). |
| `--out PATH` | Output path prefix. For each dataset writes `{base}_{DS}_table.csv` (standard classification table) and `{base}_{DS}_tests.csv` (test results) beside PATH, e.g. `--out results/run.csv` → `results/run_GSE144735_table.csv`. |
| `--no-distro` | Skip the per-feature count-distribution scan. REF and MUX tests still run, and DIST runs with minimal placeholder statistics. Much faster on large files (backed HDF5 reads). |
| `--chunk-size N` | Cells per chunk for the distribution scan (default: 500). |
| `--verbose` | Print full per-feature tables (classification, standard table, tests) to stdout. |
| `--brief` | One-line-per-finding digest (default output level). |
| `--strict` | Exit code 2 when any test is WARN (0/1 otherwise as usual). |
| `--refresh-hgnc` | Opt-in network access: query the HGNC REST API for symbols missing from the local HGNC cache before running. |

Dataset mode rules:
- **Direct file mode**: `scverify path/to/DS.h5ad` — verifies that one file in
  place. If a sibling `raw/` directory with a feature reference exists, it is
  used for MUX/ref reconciliation.
- **Base-dir mode**: `scverify --all --base-dir DIR` — every `DIR/<DS>/` with
  an `.h5ad` inside is verified, and `DIR/<DS>/raw/` is used when present.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | No FAIL; at most WARNs (unless `--strict`) |
| 1 | At least one test FAILed (curator review required; nothing was modified) |
| 2 | `--strict` and at least one WARN |

## Cache management

The `cache/` directory ships prebuilt and the tool runs offline by default:

| File | Provenance |
|------|-----------|
| `gencode_v44_*.json` (5 files) | GENCODE v44 GTF (`cache/gencode.v44.annotation.gtf.gz`) |
| `ncbi_gene_mapping.json` | NCBI `Homo_sapiens.gene_info.gz` |
| `hgnc_mapping.json` | HGNC REST snapshot |
| `adt_reference_panel.json` | GSE262381 `obsm['adt']` (137 × 26,082) |

Rebuild/refresh via `tools/fetch_caches.py` (options: `--dest DIR`,
`--gtf URL|path`, `--gene-info URL|path`, `--only gencode|ncbi|hgnc`,
`--hgnc-symbols-file FILE`):

```bash
# Rebuild GENCODE caches from the bundled GTF (offline):
python tools/fetch_caches.py --only gencode --gtf cache/gencode.v44.annotation.gtf.gz

# Rebuild the NCBI alias map from the bundled gene_info.gz (offline):
python tools/fetch_caches.py --only ncbi --gene-info cache/Homo_sapiens.gene_info.gz

# Default run: GENCODE + NCBI. If the raw files already exist in the cache
# dir, they are reused, so the default run is also fully offline.
python tools/fetch_caches.py

# Refresh the HGNC cache for a list of symbols (opt-in network step):
python tools/fetch_caches.py --only hgnc --hgnc-symbols-file symbols.txt

# Rebuild the ADT reference panel from a CITE-seq dataset:
python tools/build_adt_panel.py /path/to/cite_seq.h5ad --obsm-key adt \
    --provenance "GSE262381 obsm[adt]" \
    --out cache/adt_reference_panel.json
```

`--refresh-hgnc` on the verification CLI is the only network path inside a
normal run, and it only fills symbols missing from `hgnc_mapping.json`.

**NCBI re-parse caveat**: the alias map is a point-in-time snapshot. In the
current NCBI gene_info snapshot the mitochondrial rows (symbols `COX1`,
`COX2`, `COX3`, `ND1`) follow the nuclear rows that list the same strings as
synonyms, so a fresh last-wins re-parse resolves them to the mitochondrial
genes, while the shipped canonical cache resolves them to the nuclear genes
(e.g. `COX1`→`PTGS1`). All validated runs used the shipped canonical cache;
treat the 23 mitochondrial rRNA/tRNA keys that appear only on re-parse as
snapshot additions, not corrections.

## Worked examples

### 1. FAIL: antibody counts inside the RNA matrix

Synthetic dataset: 200 cells × 52 genes, one ADT-spelled var_name with
nonzero counts in X, an `HTO-1` obsm/var collision.

```bash
python -m scverify tmp/synth_adt_hto.h5ad --no-distro --brief
```

```
  OVERALL: FAIL (antibody counts in RNA matrix) — issue found, curator review required, nothing was removed
  ADT: ADTs are present and are placed in X/var (1 ADT-like var_names, 1 nonzero, 0 zero-count placeholders); no antibody matrix in obsm.
  MUX: multiplexing tags discovered (0 ref Multiplexing Capture, 0 HTO-ID, 1 spelling, 2 bookkeeping): 0 nonzero in X/var, 1 zero-count in X/var, 2 in obsm.
  [FAIL] antibody counts in RNA matrix — 1/1 ADT-like var_names have NONZERO counts in X — possible ADT/GEX mixing; curator review required, nothing removed by this script
  [WARN] gene-name / barcode conflicts — gene-name/feature-type conflicts (zero-count:[] unchecked:[]; var/obsm collisions:['HTO-1']); retained
```

Exit code 1.

### 2. FAIL: HTO counts inside the RNA matrix

NK2: 5569 cells × 22225 genes, 8 custom-hashtag columns (patient/sample IDs
like `NK825`) carrying nonzero RNA-matrix counts.

```bash
python -m scverify NK2.h5ad --brief
```

```
  OVERALL: FAIL (hashtag placement) — issue found, curator review required, nothing was removed
  ADT: No antibody signals detected: no ADT matrix in obsm and no ADT-like var_names in X/var.
  MUX: multiplexing tags discovered (0 ref Multiplexing Capture, 8 HTO-ID, 0 spelling, 6 bookkeeping): 8 nonzero in X/var, 0 zero-count in X/var, 0 in obsm. unsupervised discovery independently points at 6: NK1024, NK1143, NK569, NK1147, NK945, NK825.
  [FAIL] hashtag placement — 8/8 multiplexing-tag names (ref:0 htoID:8 spelling:0 bookkeeping:6) have NONZERO counts in X — possible HTO/GEX mixing; curator review required, nothing removed
  [WARN] demultiplexing consistency — demux bookkeeping inconsistencies (...)
  [WARN] unsupervised hashtag search — unsupervised counts-only discovery points at 6 hashtag-like columns (...)
```

Exit code 1.

### 3. PASS WITH NOTES: clean GEX dataset

GSE144735: 27,414 cells × 33,694 genes, pure GEX.

```bash
python -m scverify file_database/GSE144735/GSE144735.h5ad --no-distro --brief
```

```
  OVERALL: PASS WITH NOTES (ncRNA names grounded in GENCODE)
  [WARN] ncRNA names grounded in GENCODE — ncRNA-pattern names: 1987 confirmed ncRNA biotype, 54 other GENCODE biotype, 47 unconfirmed (absent from GENCODE v44; retained)
```

Exit code 0 (2 only under `--strict`).

## Performance notes

- The per-feature distribution scan reads all of X in `--chunk-size` chunks;
  on multi-GB files it can take tens of minutes (e.g. ~36 min for the 9.2 GB
  GSE178341). Use `--no-distro` for fast REF/MUX audits and run the full
  distribution pass separately when count-level evidence is needed.
- Files are always opened `backed="r"` (read-only); scVerify never writes to
  the dataset.
