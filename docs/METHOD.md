# Method

scVerify is a **read-only** auditor for single-cell `.h5ad` files. It answers
one operational question — *are antibody/hashtag counts mixed into the RNA
matrix, and do feature names resolve against authoritative references?* —
without modifying anything. Every feature is retained; findings are
classification and test results, never filters.

## Test domains

| Domain | Question | Tests |
|--------|----------|-------|
| REF  | Do feature names resolve to real genes? | T-REF-01 … T-REF-10 |
| DIST | Do count distributions look like RNA, or like antibody/hashtag counts? | T-DIST-01 … T-DIST-04 |
| MUX  | Where are the multiplexing (HTO) tags, and is the demux bookkeeping consistent? | T-MUX-01 … T-MUX-03 |

Results are `PASS`, `WARN`, or `FAIL` with a stable test ID, human-readable
check title, affected count, and up to 5 examples. Overall verdict:
`FAIL` if any test fails, else `PASS WITH NOTES` if any warn, else `CLEAN`.

## Evidence tiers (first match wins)

Feature names are resolved in strict order; the first tier that matches
supplies the resolved symbol:

1. **E1 — `E1_ensg_match`**: var-level Ensembl ID (from a var column such as
   `ensembl_id`/`gene_ids`, or the var_name itself), version suffix stripped,
   matched against the GENCODE v44 ENSG set.
2. **E2 — `E2_gencode_symbol`**: var_name is a GENCODE v44 symbol.
3. **E3 — HGNC**: resolved via the local HGNC cache.
   `E3a_hgnc_renamed` (name maps to a differing current symbol) or
   `E3b_hgnc_current` (name is already the current symbol). Cache entries
   whose value is unknown fall through to E4/E5. Offline by default;
   `--refresh-hgnc` queries the HGNC REST API for missing symbols (opt-in).
4. **E4 — `E4_ncbi_alias`**: resolved via the local NCBI Gene cache
   (alias → current symbol, from `Homo_sapiens.gene_info.gz`).
5. **E5 — `E5_unresolved_retained`**: not found in any authority. Retained,
   annotated with descriptive spelling tags only.

Biotype is looked up per tier through the GENCODE maps. A feature is
**confirmed** when it has a biotype other than `artifact`.

GENCODE v44 census used to define the ncRNA biotype set: protein_coding 20046,
lncRNA 18866, misc_RNA 2212, snRNA 1901, miRNA 1879, snoRNA 942, scaRNA 49,
rRNA 47, Mt_tRNA 22, ribozyme 8, sRNA 5, Mt_rRNA 2, vault_RNA 1, scRNA 1.
The 19-gene `artifact` biotype set is flagged by GENCODE itself (T-REF-07).

## Naming-pattern tags (spelling only)

E5 features (and, as weak fallback evidence for MUX, all features) receive
descriptive tags from a fixed regex catalogue (e.g. `adt-like`, `hto-like`,
`rfam-like`, `snorna-like`, `mir-like`, `lincrna-like`, `clone-like`,
`accession-like`, `versioned`). **Tags describe spelling only — never
biology.** "clone-like" means *the string looks like a clone-derived
identifier*; it does not mean the feature is an artifact or safe to drop. The
catalogue was built from the union of var_names across GSE132465, GSE144735,
GSE178341, GSE221575, GSE262381 (71,148 unique names) plus HGNC cache keys.

## REF tests

| ID | Check | FAIL when |
|----|-------|-----------|
| T-REF-01 | feature names are unique | duplicates exist |
| T-REF-02 | feature names are non-blank | blank names exist |
| T-REF-03 | Ensembl IDs are well-formed | ENSG-like names are malformed |
| T-REF-04 | Ensembl IDs exist in GENCODE v44 | (WARN) valid-format IDs absent from v44 |
| T-REF-05 | symbols agree with Ensembl IDs | (WARN) symbol/ENSG pairs disagree |
| T-REF-06 | ncRNA names grounded in GENCODE | (WARN) ncRNA-pattern names not confirmed by biotype |
| T-REF-07 | GENCODE-flagged artifacts scan | (WARN) features resolve to `artifact` biotype |
| T-REF-08 | reference coverage guard | (FAIL) authority caches missing/corrupt |
| T-REF-09 | feature barcodes reconciled | (FAIL) reference rows unaccounted for |
| T-REF-10 | gene-name / barcode conflicts | (FAIL) a CITE/HTO tag shares a gene name with nonzero counts; (WARN) zero-count or unchecked collisions |

## DIST tests — the empirical ADT envelope

The reference panel `cache/adt_reference_panel.json` is built from **real**
CITE-seq counts (GSE262381 `obsm['adt']`, 137 antibody columns × 26,082
cells). For each antibody it stores mean, var, overdispersion (var/mean),
frac_expressed, median, max, and the bimodality coefficient
`bc = (skew² + 1) / kurtosis` (Fisher kurtosis off); `envelope_p5_p95` holds
the p5/p50/p95 of each statistic across the 137 antibodies. Thresholds are
calibrated on those counts, not invented.

| ID | Check | Semantics |
|----|-------|-----------|
| T-DIST-01 | protein-envelope screen | Features (non-confirmed) whose mean/overdisp/frac/bc all fall inside the ADT envelope are suspects. (WARN) |
| T-DIST-02 | antibody counts in RNA matrix | ADT-spelled var_names with **nonzero** counts in X = mixing evidence. (FAIL) |
| T-DIST-03 | RNA-vs-protein sparsity check | (WARN) if the dataset's global sparsity profile is closer to the ADT panel's than RNA-like |
| T-DIST-04 | empty columns census | Informational count of zero-count columns. (PASS/WARN) |

Confirmed transcripts (E1/E2) can still be reported by T-DIST-01 for
awareness but are never counted as HTO/ADT *suspects*.

## MUX tests

Hashtag discovery is tiered: (1) feature reference file (authoritative:
`Multiplexing Capture` rows, HTO-ID), (2) var/obs/uns **bookkeeping** traces
(`hash_classification`, `htos`, `filtered_hashes`, etc.), (3) `hto-like`
spelling as a weak, explicitly-labeled fallback. Custom hashtag names
(patient/sample IDs) match none of the spelling patterns and are discoverable
only via the feature reference or bookkeeping.

| ID | Check | Semantics |
|----|-------|-----------|
| T-MUX-01 | hashtag placement | Discovered tags with **nonzero** counts in X/var → (FAIL); zero-count placeholders or clean obsm placement → PASS |
| T-MUX-02 | demultiplexing consistency | (WARN) when demux bookkeeping columns/uns lists are internally inconsistent |
| T-MUX-03 | unsupervised hashtag search | Counts-only, no names: bimodal (`bc`-based) columns that partition cells. (WARN) by design — a gene can be bimodal too |

## Residual limit (by design)

A single gene-named tag with **no** panel, ID, bookkeeping, or type trace is
indistinguishable from a genuinely bimodal gene by spelling and distribution
alone. scVerify reports every such case as a WARN with the full evidence and
leaves the call to the curator; it never drops data.

## Authority caches (offline-first)

| File | Source | Contents |
|------|--------|----------|
| `gencode_v44_symbols.json` | GENCODE v44 GTF | gene symbols |
| `gencode_v44_ensg.json` | GENCODE v44 GTF | version-stripped ENSG IDs |
| `gencode_v44_ensg_to_symbol.json` | GENCODE v44 GTF | ENSG → symbol |
| `gencode_v44_ensg_to_biotype.json` | GENCODE v44 GTF | ENSG → gene_type |
| `gencode_v44_symbol_to_biotypes.json` | GENCODE v44 GTF | symbol → biotype set |
| `hgnc_mapping.json` | HGNC REST (cached) | name → current symbol or null |
| `ncbi_gene_mapping.json` | NCBI `Homo_sapiens.gene_info.gz` | symbol + every synonym → current symbol |
| `adt_reference_panel.json` | GSE262381 `obsm['adt']` | empirical ADT envelope |
| `gencode.v44.annotation.gtf.gz` | EBI | raw source for offline re-parse |
| `Homo_sapiens.gene_info.gz` | NCBI | raw source for offline re-parse |

Network access is required **only** for explicit refreshes (`--refresh-hgnc`,
`tools/fetch_caches.py` downloads). Re-parsing note: the NCBI alias map is
snapshot-dependent; in current gene_info snapshots the mitochondrial gene
rows (symbols `COX1`, `COX2`, `COX3`, `ND1`, …) appear after the nuclear
alias rows that list them as synonyms, so a last-wins re-parse resolves those
symbols to the mitochondrial genes. The shipped canonical cache predates that
snapshot and resolves them to the nuclear genes (e.g. `COX1`→`PTGS1`), which
is the behavior all validated runs used.
