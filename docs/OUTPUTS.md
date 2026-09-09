# Outputs

All outputs are **additive reports**. scVerify never writes to the dataset
file itself (opened `backed="r"`).

## Standard output (stdout)

```
GENCODE v44: <n_symbols> symbols, <n_ensg> ENSG IDs, maps: ...
Caches: HGNC=<n> entries, NCBI=<n> entries (offline; use --refresh-hgnc to query HGNC)
ADT panel: <provenance>
How to read: per dataset, OVERALL says CLEAN / PASS WITH NOTES / FAIL. ...
Feature reference for <DS>: <n> rows (<provenance>)   # or "none (...)"

======================================================================
  <DS>
======================================================================
  <n_cells> cells x <n_genes> genes (read-only, backed)
  var feature_types: {...}                    # only when var carries types
  OVERALL: <CLEAN | PASS WITH NOTES | FAIL> (<failing/warn check titles>)
  ADT: <placement sentence>
  MUX: <discovery sentence with source breakdown + unsupervised pointers>
  [WARN|FAIL] <check title> — <details, counts, up to 5 examples>
  ...

======================================================================
  SUMMARY
======================================================================
  <DS>: <n_cells> cells x <n_genes> genes | <FAIL:checks | WARN:checks | ALL PASS>
    OVERALL: <verdict>
    <ADT placement sentence>
    <MUX placement sentence>
```

With `--brief` (default), only WARN/FAIL checks print details. With
`--verbose`, full per-feature tables (classification, standard table, tests)
are printed as pandas DataFrames.

`OVERALL` per dataset: `FAIL (<check>)` if any test failed, else
`PASS WITH NOTES (<check>)` if any warned, else `CLEAN`.

## CSV outputs (`--out`)

`--out PATH` writes, per dataset, two CSVs named from the PATH prefix
(`base` = PATH without extension):

### `<base>_<DS>_table.csv` — standard classification table

| Column | Meaning |
|--------|---------|
| `category` | Human-readable label: one row per evidence tier (`match · Ensembl ID`, `match · GENCODE symbol`, `match · HGNC alias (renamed)`, `match · HGNC current symbol`, `match · NCBI alias`, `unresolved · kept`), then `unresolved · <tag> spelling` rows (E5 tag breakdown), `biotype: <type>` rows (E1/E2 breakdown), and `reference: <type>` rows (feature-reference types) |
| `basis` | The authority/source that established the row |
| `count` | Number of features in the row |
| `pct` | Percentage of all features (2 dp) |
| `examples` | Up to 5 example names, `; `-joined. HGNC-renamed rows use `<name>-><resolved>` form |
| `disposition` | Always a retain/report note — never a drop |

### `<base>_<DS>_tests.csv` — test results

| Column | Meaning |
|--------|---------|
| `test` | Stable test ID (`T-REF-01` … `T-MUX-03`) |
| `check` | Human-readable check title |
| `domain` | `REF`, `DIST`, or `MUX` |
| `result` | `PASS`, `WARN`, or `FAIL` |
| `summary` | One-line result description with counts |
| `n_affected` | Number of features the finding covers |
| `examples` | Up to 5 example names, `; `-joined |

## Per-feature classification (internal; printed with `--verbose`)

One row per feature, columns:

`var_name`, `evidence` (tier code), `resolved_symbol` (current symbol from
HGNC/NCBI; the name itself for E1/E2; None for E5), `ensg`, `biotype`,
`ref_feature_type` (authoritative type from the feature reference file, or
None), `primary_tag`, `all_tags` (`|`-joined), `confirmed` (True when an
authority identifies a known gene transcript with a non-`artifact` biotype).

## Cache files

| File | Format |
|------|--------|
| `gencode_v44_symbols.json` | JSON array of gene symbols |
| `gencode_v44_ensg.json` | JSON array of version-stripped ENSG IDs |
| `gencode_v44_ensg_to_symbol.json` | JSON object ENSG → symbol |
| `gencode_v44_ensg_to_biotype.json` | JSON object ENSG → gene_type |
| `gencode_v44_symbol_to_biotypes.json` | JSON object symbol → [biotypes] |
| `hgnc_mapping.json` | JSON object name → current symbol (or null) |
| `ncbi_gene_mapping.json` | JSON object symbol/synonym → current symbol |
| `adt_reference_panel.json` | JSON: `provenance`, `n_cells`, `n_antibodies`, `features` (per-antibody `mean`, `var`, `overdisp`, `frac_expressed`, `bc`, `median`, `max`), `envelope_p5_p95` (p5/p50/p95 of each statistic) |
| `{DS}_feature_reference.json` | Per-dataset cache of the parsed feature reference (feature_id → feature_type), built on demand from the dataset's `raw/` directory |

## Exit codes

`0` clean or WARN-only; `1` at least one FAIL; `2` WARN-only with `--strict`.
