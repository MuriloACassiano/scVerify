"""Verification std codes: ADT-in-GEX checks + feature-name reference tests.

READ-ONLY. Nothing here modifies, filters, or deletes data. Functions classify,
test, and report only.

Test domains:
  REF  — names/symbols/ncRNA identifiers checked against GENCODE v44 (symbols,
         ENSG IDs, gene_type biotypes incl. a 19-gene "artifact" set), the local
         HGNC cache, and the local NCBI Gene cache. Includes a symbol<->ENSG
         cross-column consistency check wherever var carries both, plus
         feature-barcode reconciliation against the feature reference file.
  DIST — per-feature count distributions tested against an empirical ADT
         reference panel built from REAL CITE-seq counts (GSE262381 obsm['adt'],
         137 antibodies; cache/adt_reference_panel.json). Envelope and
         thresholds are calibrated on those counts, not invented.
  MUX  — multiplexing (HTO) tag discovery and placement. Custom hashtag names
         (patient/sample IDs) are discovered ONLY via discover_feature_reference
         (raw/ feature_ref.csv / features.tsv / 10x H5); hto-like spelling is
         a weak fallback that is reported as such.

Evidence tiers (authoritative, applied in order — first match wins):
  E1  var-level Ensembl ID  — ENSG ID taken from a var column (ensembl_id/gene_ids)
                              or from var_name itself, matched against the
                              GENCODE v44 ENSG set (version suffix stripped).
  E2  GENCODE symbol        — var_name matched against GENCODE v44 symbols.
  E3  HGNC alias resolution — var_name resolved via the local HGNC cache
                              (cache/hgnc_mapping.json): renamed / same / unknown.
                              Offline by default; network refresh only via
                              refresh_hgnc_cache (opt-in).
  E4  NCBI alias resolution — var_name resolved via the local NCBI Gene cache
                              (cache/ncbi_gene_mapping.json), built from
                              Homo_sapiens.gene_info.gz. Offline.
  E5  Unresolved, retained  — not found in any authority above. Kept as-is with
                              descriptive naming-pattern tags (see PATTERNS).
                              Tags describe SPELLING ONLY, never biology.

Naming-pattern tags are descriptive, never biological classes. A tag such as
"clone-like" means "the string looks like a clone-derived identifier". It does
NOT mean the feature is an artifact, a contaminant, or safe to drop.

CLI entry point: `python -m scverify` (or the `scverify` console script after
`pip install -e .`).
"""

import json
import os
import re

import numpy as np
import pandas as pd

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cache")
HGNC_API = "https://rest.genenames.org"

# GENCODE v44 gene_type values treated as non-coding RNA for the ncRNA test.
# (Auditable: full gene_type census was taken from cache/gencode.v44.annotation.gtf.gz;
#  counts: protein_coding 20046, lncRNA 18866, misc_RNA 2212, snRNA 1901,
#  miRNA 1879, snoRNA 942, scaRNA 49, rRNA 47, Mt_tRNA 22, ribozyme 8,
#  sRNA 5, Mt_rRNA 2, vault_RNA 1, scRNA 1.)
NCRNA_BIOTYPES = frozenset({
    "lncRNA", "miRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA", "rRNA",
    "Mt_tRNA", "Mt_rRNA", "ribozyme", "sRNA", "vault_RNA", "scRNA",
})

# var_name spellings that motivate an ncRNA-authority check (descriptive trigger only).
NCRNA_PATTERN_TAGS = frozenset({"snorna-like", "mir-like", "lincrna-like", "rfam-like"})

# ============================================================================
# Comprehensive naming-pattern catalogue (descriptive ONLY — spelling, not biology)
# Each entry: (tag, compiled regex, note). Order matters: first match wins for
# the *primary* tag; all matches are also listed in `all_tags`.
# This list was built from the union of var_names across GSE132465, GSE144735,
# GSE178341, GSE221575 and GSE262381 (71,148 unique names) plus HGNC cache keys.
# ============================================================================

PATTERNS = [
    # --- Ensembl identifiers used as var_names (evidence tier E1 handles these) ---
    ("ensembl-id-like",
     re.compile(r"^ENSG\d{11}(\.\d+)?$"),
     "Spelling matches an Ensembl gene ID, with or without version suffix."),

    # --- Antibody-derived-tag spellings (CITE-seq feature-reference style) ---
    ("adt-like",
     re.compile(r"^(ADT-.+|Hu(\.|Ms).+|Isotype_.+|.*-Isotype.*)$"),
     "Spelling resembles an antibody-derived tag (e.g. ADT-CD4, Hu.CD8, "
     "Isotype_*). Says nothing about counts or storage; check both."),

    # --- Multiplexing/hashtag spellings (CITE-seq/CellPlex feature-reference style) ---
    ("hto-like",
     re.compile(r"^(HTO([-_.\s].*)?$|[Hh]ashtag|CMO([-_.].*)?$|MULTI[-_. ]?seq"
                r"|Lipid[-_. ]?tag|CellPlex|hto_"
                r"|.*TotalSeq.*HTO.*|.*[Hh]ashtag.*|.*HTO\d+[A-Z]*$)"),
     "Spelling resembles a multiplexing/hashtag tag (e.g. HTO-1, Hashtag-A, "
     "CMO-..., TotalSeqB0251_HTO1). Custom hashtag names (patient/sample IDs "
     "such as NK825) do NOT match this — those are discoverable only via the "
     "feature reference file or var/obs/uns bookkeeping. Spelling only; "
     "check counts and source."),

    # --- Rfam / small-RNA family identifiers ---
    ("rfam-like",
     re.compile(r"^RF\d+(_\d+)?$"),
     "Spelling matches an Rfam family ID (e.g. RF00001, RF00003_11)."),

    # --- Small non-coding RNA spellings ---
    ("snorna-like",
     re.compile(r"^(snoU\d+(-\d+)?|SNORD\d+.*|SNORA\d+.*|RNU\d+.*|RN7[A-Z]*\d*P?\d*"
                r"|Y_RNA.*|U\d+-\d+|5S_rRNA|5_8S_rRNA|7SK)$"),
     "Spelling resembles a small ncRNA (snoRNA, snRNA, Y RNA, U3 variants)."),

    # --- MicroRNA / lncRNA spellings (many are also GENCODE symbols; E2 wins) ---
    ("mir-like",
     re.compile(r"^MIR\d+.*$"),
     "Spelling resembles a microRNA host or precursor symbol."),
    ("lincrna-like",
     re.compile(r"^(LINC\d+.*|.*-AS\d*$|.*-HG$|HOTAIR|MALAT1|NEAT1|XIST|FAM138A)$"),
     "Spelling resembles a long non-coding RNA symbol."),

    # --- Clone-derived identifiers (BAC/PAC/fosmid/cosmid libraries) ---
    ("clone-like",
     re.compile(r"^((RP\d*|CT[ABCD]|CTC|XXbac|XXyac|LLNLR|KB|CH\d+|LA\d+c|GS\d+)"
                r"-[A-Za-z0-9]+\.\d+|XXbac-[A-Za-z0-9]+\.\d+)$"),
     "Spelling matches a genomic-clone-derived identifier "
     "(e.g. RP11-543A18.1, CTD-3074O7.5, XXbac-B135H6.18)."),

    # --- GenBank/EMBL/DDBJ accession-style identifiers (versioned) ---
    ("accession-like",
     re.compile(r"^([A-Z]{1,2}\d{4,8}|[A-Z]{3,4}\d{2,8}|[A-Z]{2}\d+)"
                r"[A-Za-z0-9]*(\.\d+)?$"),
     "Spelling matches a sequence-accession style token such as AC087072.1, "
     "AF001548.2, AL627309.1, AP006222.2, BX004987.1, Z68694.1. "
     "Covers 1-4 letter prefixes + digits, optional version suffix."),

    # --- Generic version suffix (e.g. SYMBOL.1, ORF.2) ---
    ("versioned",
     re.compile(r"^.+\.\d+$"),
     "Name carries a dot-number version suffix (e.g. AL627309.1, AKR1C3.1). "
     "The suffix is part of the annotation string, not a separate gene."),

    # --- Double-underscore / underscore-digit suffixes ---
    ("underscore-suffixed",
     re.compile(r"^.+__.*$|^.+_\d+$"),
     "Name carries __ or _N suffix (e.g. RP11-544L8__B.4, RF00001_2)."),

    # --- Hyphenated symbols (covers readthroughs, antisense, pseudogenes) ---
    ("hyphenated",
     re.compile(r"^.+-.+$"),
     "Name contains a hyphen (e.g. A1BG-AS1, MT-ND1, HLA-DRA)."),

    # --- Fallback: plain symbol ---
    ("plain-symbol",
     re.compile(r"^.+$"),
     "None of the above spellings; treated as a plain symbol string."),
]


def tag_name(name):
    """Return (primary_tag, all_tags) for a var_name. Descriptive only."""
    hits = [tag for tag, rx, _ in PATTERNS if rx.match(name)]
    primary = hits[0] if hits else "plain-symbol"
    return primary, hits


# ============================================================================
# Authority caches (offline-first)
# ============================================================================

def _load_json(path):
    with open(path) as f:
        return json.load(f)


def load_gencode_symbols(cache_dir=None):
    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache_file = os.path.join(cache_dir, "gencode_v44_symbols.json")
    if not os.path.exists(cache_file):
        raise FileNotFoundError(
            f"GENCODE symbol cache not found at {cache_file}."
        )
    data = _load_json(cache_file)
    return set(data)


def load_gencode_ensg(cache_dir=None):
    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache_file = os.path.join(cache_dir, "gencode_v44_ensg.json")
    if not os.path.exists(cache_file):
        return set()
    return set(_load_json(cache_file))


def load_gencode_maps(cache_dir=None):
    """Load ENSG->symbol, ENSG->biotype, symbol->biotypes maps (empty dicts if absent)."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    out = {}
    for key, fn in [("ensg_to_symbol", "gencode_v44_ensg_to_symbol.json"),
                    ("ensg_to_biotype", "gencode_v44_ensg_to_biotype.json"),
                    ("symbol_to_biotypes", "gencode_v44_symbol_to_biotypes.json")]:
        p = os.path.join(cache_dir, fn)
        out[key] = _load_json(p) if os.path.exists(p) else {}
    return out


def load_hgnc_cache(cache_dir=None):
    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache_file = os.path.join(cache_dir, "hgnc_mapping.json")
    if os.path.exists(cache_file):
        return _load_json(cache_file)
    return {}


def load_ncbi_cache(cache_dir=None):
    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache_file = os.path.join(cache_dir, "ncbi_gene_mapping.json")
    if os.path.exists(cache_file):
        return _load_json(cache_file)
    return {}


def load_adt_panel(cache_dir=None):
    """Load the empirical ADT reference panel (None if absent)."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    p = os.path.join(cache_dir, "adt_reference_panel.json")
    return _load_json(p) if os.path.exists(p) else None


def discover_feature_reference(ds, base_dir, cache_dir=None, refresh=False,
                               raw_dir=None):
    """Authoritative feature-type discovery from raw/ files (read-only).

    Custom hashtag/antibody names (patient/sample IDs, TotalSeq codes) cannot
    be discovered by spelling — they are discovered HERE, from the files the
    counts were built with. Searches file_database/<DS>/raw/ for (in order):
      1. *feature_ref.csv — Cell Ranger feature reference
         (id,name,...,feature_type). Covers custom names automatically.
      2. *features.tsv[.gz] — matrix feature tables (id, name, feature_type).
         Old *genes.tsv files carry no type column and are reported as such.
      3. *.h5 — 10x H5 matrices (matrix/features/{name,feature_type}).
      4. *.tar — archives are listed and matching members (*feature_ref.csv,
         *features.tsv*, *genes.tsv*) are parsed in-memory via tarfile
         (no disk writes).
    Result {feature_name: feature_type} is cached to
    cache/<DS>_feature_reference.json with provenance. Returns (ref_map or
    None, provenance string). Files outside cache/ are never written.
    """
    import csv
    import fnmatch
    import gzip
    import tarfile

    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache_file = os.path.join(cache_dir, f"{ds}_feature_reference.json")
    if os.path.exists(cache_file) and not refresh:
        cached = _load_json(cache_file)
        return cached.get("features") or None, cached.get("provenance", "")

    raw_dir = os.path.join(base_dir, ds, "raw") if raw_dir is None else raw_dir
    seen, ref_map, provenance = [], None, ""
    if not raw_dir or not os.path.isdir(raw_dir):
        provenance = "no raw/ directory"
    else:
        entries = sorted(os.listdir(raw_dir))
        seen = entries

        def _adopt(mapping, prov):
            nonlocal ref_map, provenance
            if mapping and ref_map is None:
                ref_map, provenance = mapping, prov

        # 1. Cell Ranger feature reference CSVs
        for fn in entries:
            if fnmatch.fnmatch(fn, "*feature_ref.csv"):
                try:
                    with open(os.path.join(raw_dir, fn), newline="") as fh:
                        rows = list(csv.DictReader(fh))
                    if rows and "name" in rows[0] and "feature_type" in rows[0]:
                        _adopt({r["name"]: r["feature_type"] for r in rows
                                if r.get("name")},
                               f"{fn} (Cell Ranger feature reference)")
                except (OSError, csv.Error):
                    continue

        # 2. Matrix feature tables (plain or gzipped)
        for fn in entries:
            if fnmatch.fnmatch(fn, "*features.tsv*") and ref_map is None:
                p = os.path.join(raw_dir, fn)
                try:
                    opener = gzip.open if fn.endswith(".gz") else open
                    mapping, has_type = {}, False
                    with opener(p, "rt") as fh:
                        for line in fh:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) >= 3 and parts[2].strip():
                                mapping[parts[1]] = parts[2]
                                has_type = True
                    if has_type:
                        _adopt(mapping, f"{fn} (matrix features.tsv)")
                except OSError:
                    continue

        # 3. 10x H5 matrices
        for fn in entries:
            if fnmatch.fnmatch(fn, "*.h5") and ref_map is None:
                try:
                    import h5py
                    with h5py.File(os.path.join(raw_dir, fn), "r") as h5:
                        grp = h5.get("matrix/features", h5.get("features", None))
                        if grp is None or "name" not in grp \
                                or "feature_type" not in grp:
                            continue
                        names = [n.decode() if isinstance(n, bytes) else str(n)
                                 for n in grp["name"][:]]
                        ftypes = [t.decode() if isinstance(t, bytes) else str(t)
                                  for t in grp["feature_type"][:]]
                    _adopt(dict(zip(names, ftypes)),
                           f"{fn} (10x H5 matrix/features)")
                except (OSError, ImportError, KeyError):
                    continue

        # 4. Look inside tar archives (in-memory; no disk writes)
        for fn in entries:
            if fnmatch.fnmatch(fn, "*.tar") and ref_map is None:
                try:
                    with tarfile.open(os.path.join(raw_dir, fn), "r") as tf:
                        members = [m for m in tf.getmembers()
                                   if m.isfile() and (
                                       fnmatch.fnmatch(m.name, "*feature_ref.csv")
                                       or fnmatch.fnmatch(m.name, "*features.tsv*"))]
                        for m in members:
                            fh = tf.extractfile(m)
                            if fh is None:
                                continue
                            if m.name.endswith("feature_ref.csv"):
                                text = fh.read().decode()
                                rows = list(csv.DictReader(text.splitlines()))
                                if rows and "name" in rows[0] \
                                        and "feature_type" in rows[0]:
                                    _adopt(
                                        {r["name"]: r["feature_type"] for r in rows
                                         if r.get("name")},
                                        f"{fn}:{m.name} (feature reference in TAR)")
                                    break
                            else:  # features.tsv inside tar
                                content = fh.read().decode()
                                mapping = {}
                                for line in content.splitlines():
                                    parts = line.split("\t")
                                    if len(parts) >= 3 and parts[2].strip():
                                        mapping[parts[1]] = parts[2]
                                if mapping:
                                    _adopt(mapping,
                                           f"{fn}:{m.name} (features.tsv in TAR)")
                                    break
                except (OSError, tarfile.TarError):
                    continue

        if ref_map is None:
            genes_only = [f for f in seen if "genes.tsv" in f]
            provenance = ("no typed feature source in raw/ "
                          f"(files seen: {', '.join(seen[:8])}"
                          f"{'...' if len(seen) > 8 else ''})"
                          + (f"; genes.tsv present (no type column): "
                             f"{', '.join(genes_only)}" if genes_only else ""))

    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_file, "w") as fh:
        json.dump({"provenance": provenance, "features": ref_map or {}}, fh)
    return ref_map, provenance


def load_all_caches(cache_dir=None):
    """Load every authority cache + the ADT panel. Returns a dict."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    return {
        "gencode_symbols": load_gencode_symbols(cache_dir),
        "gencode_ensg": load_gencode_ensg(cache_dir),
        "gmaps": load_gencode_maps(cache_dir),
        "hgnc": load_hgnc_cache(cache_dir),
        "ncbi": load_ncbi_cache(cache_dir),
        "adt_panel": load_adt_panel(cache_dir),
    }


def discover_var_feature_types(var_names, var_feature_types):
    """Authoritative in-matrix feature types from a var column (read-only view).

    Returns {name: type} for non-empty entries, else {}. The var feature_types
    column describes exactly what is in THIS X — including custom hashtag names
    (e.g. NK825 -> 'Antibody Capture', grounded further by TotalSeq HTO IDs in
    a gene_ids column) with no external files needed.
    """
    if var_feature_types is None:
        return {}
    out = {}
    for n, t in zip(var_names, list(var_feature_types)):
        t = str(t).strip() if t is not None else ""
        if t and t.lower() != "nan":
            out[n] = t
    return out


BOOKKEEPING_CLASS_KEYWORDS = frozenset({
    "hash", "hto", "demux", "multiplex", "souporcell", "soupor", "vireo",
    "demuxlet", "freemuxlet", "hashtag", "lipid", "cmo", "cellplex", "multi",
})

BOOKKEEPING_LIST_KEYWORDS = frozenset({
    "hto", "hash", "multiplex", "demux", "singlet", "donor", "barcode",
})


def _parse_uns_tag_list(value):
    """Parse an uns entry into a list of tag strings, or None if not list-like."""
    import numpy as np
    import pandas as pd
    if isinstance(value, (list, tuple, np.ndarray, pd.Index)):
        items = [str(x).strip() for x in list(value)]
        items = [x for x in items if x and x.lower() != "nan"]
        return items or None
    if isinstance(value, str):
        items = re.findall(r"'([^']+)'|\"([^\"]+)\"", value)
        flat = [a or b for a, b in items]
        if flat:
            return flat
        return [x for x in re.split(r"[\s,;|]+", value.strip("[]()")) if x] or None
    return None


def summarize_mux_bookkeeping(adata):
    """Demux bookkeeping evidence from obs/uns (read-only; descriptive).

    Detects the traces a hashtag workflow leaves behind, including custom
    names: per-cell classification columns, per-tag count columns shared
    between obs and var, and uns tag lists (kept vs filtered). Returns a dict
    with keys: class_column, class_values, tag_obs_var_overlap, uns_hto_lists,
    uns_filtered_lists. Empty findings are reported as such — absence of
    bookkeeping is itself a stated limit for custom-name discovery.
    """
    obs_cols = list(adata.obs.columns)
    var_set = set(adata.var_names.tolist())

    def _is_class_col(c):
        cl = c.lower()
        return any(k in cl for k in BOOKKEEPING_CLASS_KEYWORDS)

    class_col = next((c for c in obs_cols if _is_class_col(c)), None)
    class_values = []
    if class_col is not None:
        try:
            class_values = sorted({str(x) for x in adata.obs[class_col].tolist()
                                   if str(x).strip()})
        except (KeyError, TypeError, ValueError):
            class_values = []

    tag_overlap = sorted({c for c in obs_cols if c in var_set})

    uns_hto, uns_filt = {}, {}
    try:
        uns_items = list(adata.uns.items())
    except (AttributeError, TypeError):
        uns_items = []
    for k, v in uns_items:
        kl = str(k).lower()
        if "color" in kl:
            continue
        if not any(w in kl for w in BOOKKEEPING_LIST_KEYWORDS):
            continue
        items = _parse_uns_tag_list(v)
        if not items:
            continue
        if any(w in kl for w in ("filter", "remove", "exclud", "discard",
                                 "doublet", "negative")):
            uns_filt[str(k)] = items
        else:
            uns_hto[str(k)] = items

    return {
        "class_column": class_col,
        "class_values": class_values,
        "tag_obs_var_overlap": tag_overlap,
        "uns_hto_lists": uns_hto,
        "uns_filtered_lists": uns_filt,
    }


def load_all_caches(cache_dir=None):
    """Load every authority cache + the ADT panel. Returns a dict."""
    if cache_dir is None:
        cache_dir = CACHE_DIR
    gmaps = load_gencode_maps(cache_dir)
    return {
        "gencode_symbols": load_gencode_symbols(cache_dir),
        "gencode_ensg": load_gencode_ensg(cache_dir),
        "gmaps": gmaps,
        "hgnc": load_hgnc_cache(cache_dir),
        "ncbi": load_ncbi_cache(cache_dir),
        "adt_panel": load_adt_panel(cache_dir),
    }


ENSG_RE = re.compile(r"(ENSG\d{11})")
ENSG_LIKE_RE = re.compile(r"^ENSG")
ENSG_STRICT_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$")

# Hashtag-ID patterns for feature-ID columns (gene_ids / feature id). Verified
# zero-collision against all 61,228 GENCODE v44 symbols (the single HTO-bearing
# symbol, CHTOP, carries no trailing digit and never matches).
HTO_ID_RE = re.compile(r"TotalSeq.*HTO|HTO[-_.\s]?\d|[Hh]ashtag")


def extract_ensg(text):
    """Extract a bare ENSG ID (no version) from any string, or None."""
    if text is None:
        return None
    m = ENSG_RE.search(str(text))
    return m.group(1) if m else None


def find_ensg_column(var_columns):
    for col in ["ensembl_id", "gene_ids", "gene_id", "ensembl_gene_id"]:
        if col in var_columns:
            return col
    return None


def find_id_column(var_columns):
    """First var column usable as a feature-ID source (IDs ground custom names)."""
    for col in ["gene_ids", "gene_id", "ensembl_id", "ensembl_gene_id",
                "feature_id", "id"]:
        if col in var_columns:
            return col
    return None


# ============================================================================
# Evidence-tier classification (mutually exclusive; first match wins)
# ============================================================================

EVIDENCE_LABELS = [
    "E1_ensg_match",
    "E2_gencode_symbol",
    "E3a_hgnc_renamed",
    "E3b_hgnc_current",
    "E4_ncbi_alias",
    "E5_unresolved_retained",
]

# Human-readable names. Stable IDs (T-REF-01, E1_ensg_match, ...) stay in the
# `test`/`evidence` columns for scripts; the `check`/`category` columns and
# all verdict lines use these titles.
TEST_TITLES = {
    "T-REF-01": "feature names are unique",
    "T-REF-02": "feature names are non-blank",
    "T-REF-03": "Ensembl IDs are well-formed",
    "T-REF-04": "Ensembl IDs exist in GENCODE v44",
    "T-REF-05": "symbols agree with Ensembl IDs",
    "T-REF-06": "ncRNA names grounded in GENCODE",
    "T-REF-07": "GENCODE-flagged artifacts scan",
    "T-REF-08": "reference coverage guard",
    "T-REF-09": "feature barcodes reconciled",
    "T-REF-10": "gene-name / barcode conflicts",
    "T-DIST-01": "protein-envelope screen",
    "T-DIST-02": "antibody counts in RNA matrix",
    "T-DIST-03": "RNA-vs-protein sparsity check",
    "T-DIST-04": "empty columns census",
    "T-MUX-01": "hashtag placement",
    "T-MUX-02": "demultiplexing consistency",
    "T-MUX-03": "unsupervised hashtag search",
}

EVIDENCE_TITLES = {
    "E1_ensg_match": "match · Ensembl ID",
    "E2_gencode_symbol": "match · GENCODE symbol",
    "E3a_hgnc_renamed": "match · HGNC alias (renamed)",
    "E3b_hgnc_current": "match · HGNC current symbol",
    "E4_ncbi_alias": "match · NCBI alias",
    "E5_unresolved_retained": "unresolved · kept",
}

EV_SHORT = {
    "E1_ensg_match": "Ensembl-ID",
    "E2_gencode_symbol": "GENCODE-symbol",
    "E3a_hgnc_renamed": "HGNC-alias",
    "E3b_hgnc_current": "HGNC-current",
    "E4_ncbi_alias": "NCBI-alias",
    "E5_unresolved_retained": "unresolved",
}


def ev_display(code):
    """Short human-readable evidence label for inline messages."""
    return EV_SHORT.get(code, code)


def test_title(test_id):
    """Human-readable check title for a stable test ID."""
    return TEST_TITLES.get(test_id, test_id)


def _missing(v):
    """Missing-value check (None/nan/pd.NA) stable across pandas 2 and 3.

    pandas >= 3 stores missing entries of the default string dtype as nan,
    which is truthy — a bare `if v:` test silently misclassifies them.
    """
    return v is None or bool(pd.isna(v))


def lookup_biotype(evidence, name, resolved, ensg, gmaps):
    """GENCODE gene_type(s) for a classified feature, or None if unmapped."""
    e2s = gmaps.get("ensg_to_biotype", {})
    s2b = gmaps.get("symbol_to_biotypes", {})
    if evidence == "E1_ensg_match" and ensg:
        return e2s.get(ensg)
    if evidence == "E2_gencode_symbol":
        bts = s2b.get(name)
        return "|".join(bts) if bts else None
    if evidence in ("E3a_hgnc_renamed", "E3b_hgnc_current", "E4_ncbi_alias") \
            and resolved:
        bts = s2b.get(resolved)
        return "|".join(bts) if bts else None
    return None


def classify_features(var_names, var_ensg_list, gencode_symbols, gencode_ensg,
                       hgnc_cache, ncbi_cache, gmaps, ref_map=None):
    """Classify each var_name into exactly one evidence tier.

    Returns a DataFrame with columns:
      var_name, evidence, resolved_symbol, ensg, biotype, ref_feature_type,
      primary_tag, all_tags, confirmed
    `resolved_symbol` is the current symbol from HGNC/NCBI when applicable,
    else the var_name itself (E1/E2) or None (E5). `ref_feature_type` is the
    authoritative type from the feature reference file (None if undiscovered);
    it is the only source that covers custom hashtag/antibody names.
    `confirmed` is True when an authority identifies the feature as a known
    gene transcript (GENCODE biotype present and not 'artifact', any tier).
    Confirmed transcripts are never HTO/ADT *suspects* — T-DIST-01 skips them
    and discovery annotates (never convicts) them — but corroborated
    distribution evidence is never concealed: a confirmed gene that is BOTH
    protein-distributed AND partitions with a hashtag panel is reported for
    awareness (possible gene-named tag). Authoritative conflicts (e.g.
    Antibody Capture type on a confirmed gene symbol) surface in
    T-MUX-01/T-REF-09 regardless.
    """
    rows = []
    for i, name in enumerate(var_names):
        ensg = var_ensg_list[i] if var_ensg_list else extract_ensg(name)
        primary_tag, all_tags = tag_name(name)

        if ensg and ensg in gencode_ensg:
            ev, res = "E1_ensg_match", name
        elif name in gencode_symbols:
            ev, res = "E2_gencode_symbol", name
        elif name in hgnc_cache and hgnc_cache[name] not in (None, name):
            ev, res = "E3a_hgnc_renamed", hgnc_cache[name]
        elif name in hgnc_cache and hgnc_cache[name] == name:
            ev, res = "E3b_hgnc_current", name
        elif name in ncbi_cache:
            ev, res = "E4_ncbi_alias", ncbi_cache[name]
        else:
            ev, res = "E5_unresolved_retained", None

        bt = lookup_biotype(ev, name, res, ensg, gmaps)
        rows.append((name, ev, res, ensg, bt,
                     ref_map.get(name) if ref_map else None,
                     primary_tag, "|".join(all_tags),
                     bool(bt and bt != "artifact")))

    return pd.DataFrame(rows, columns=["var_name", "evidence",
                                       "resolved_symbol", "ensg", "biotype",
                                       "ref_feature_type",
                                       "primary_tag", "all_tags", "confirmed"])


def build_standard_table(class_df, ref_provenance=""):
    """Build the standard per-dataset classification table.

    One row per evidence tier + tag breakdown of E5 + biotype breakdown of
    resolved (E1+E2) features + feature-reference breakdown (authoritative,
    covers custom names). Columns:
      category, basis, count, pct, examples, disposition
    Disposition is always retain/report — never drop.
    """
    total = len(class_df)
    basis = {
        "E1_ensg_match": "ENSG ID (var column or var_name) in GENCODE v44 ENSG set",
        "E2_gencode_symbol": "var_name in GENCODE v44 symbol set",
        "E3a_hgnc_renamed": "HGNC cache: alias -> differing current symbol",
        "E3b_hgnc_current": "HGNC cache: symbol is current (same name)",
        "E4_ncbi_alias": "NCBI Gene cache: alias -> symbol",
        "E5_unresolved_retained": "No authority match; retained with pattern tags",
    }
    out_rows = []
    for ev in EVIDENCE_LABELS:
        sub = class_df[class_df["evidence"] == ev]
        n = len(sub)
        if ev.startswith("E3a"):
            ex = [f"{r.var_name}->{r.resolved_symbol}"
                  for r in sub.head(5).itertuples()]
        else:
            ex = sub["var_name"].head(5).tolist()
        out_rows.append({
            "category": EVIDENCE_TITLES[ev],
            "basis": basis[ev],
            "count": n,
            "pct": round(100.0 * n / max(total, 1), 2),
            "examples": "; ".join(ex),
            "disposition": "retain; no filtering applied",
        })

    # Tag breakdown of E5 only (descriptive, not a re-classification)
    e5 = class_df[class_df["evidence"] == "E5_unresolved_retained"]
    if len(e5):
        for tag, n in e5["primary_tag"].value_counts().items():
            ex = e5[e5["primary_tag"] == tag]["var_name"].head(5).tolist()
            out_rows.append({
                "category": f"unresolved · {tag} spelling",
                "basis": "naming-pattern tag within unresolved rows (spelling only, not a class)",
                "count": int(n),
                "pct": round(100.0 * n / max(total, 1), 2),
                "examples": "; ".join(ex),
                "disposition": "retain; tag is descriptive only",
            })

    # Biotype breakdown of authority-resolved (E1+E2) features
    resolved = class_df[class_df["evidence"].isin(["E1_ensg_match",
                                                   "E2_gencode_symbol"])]
    if len(resolved):
        bio = resolved["biotype"].fillna("biotype_unknown").value_counts()
        for bt, n in bio.items():
            ex = resolved[resolved["biotype"].fillna("biotype_unknown") == bt][
                "var_name"].head(5).tolist()
            out_rows.append({
                "category": f"biotype: {bt}",
                "basis": "GENCODE v44 gene_type for E1/E2 features",
                "count": int(n),
                "pct": round(100.0 * n / max(total, 1), 2),
                "examples": "; ".join(ex),
                "disposition": "retain; biotype is annotation, not a filter",
            })

    # Feature-reference breakdown (authoritative; covers custom names)
    if "ref_feature_type" in class_df.columns:
        refd = class_df[class_df["ref_feature_type"].notna()]
        if len(refd):
            for ft, n in refd["ref_feature_type"].value_counts().items():
                ex = refd[refd["ref_feature_type"] == ft]["var_name"].head(5).tolist()
                out_rows.append({
                    "category": f"reference: {ft}",
                    "basis": f"feature reference ({ref_provenance or 'unknown provenance'})",
                    "count": int(n),
                    "pct": round(100.0 * n / max(total, 1), 2),
                    "examples": "; ".join(ex),
                    "disposition": "retain; reference type is annotation, not a filter",
                })

    return pd.DataFrame(out_rows, columns=["category", "basis", "count", "pct",
                                           "examples", "disposition"])
# ============================================================================
# Multiplexing (HTO) tests + feature-barcode reconciliation
# Custom hashtag names (patient/sample IDs) are discovered ONLY via the
# feature reference file (discover_feature_reference); spelling (hto-like)
# is a weak fallback. Results: PASS / WARN / FAIL, same semantics as REF/DIST.
# ============================================================================

NON_GEX_FEATURE_TYPES = frozenset({
    "Antibody Capture", "Multiplexing Capture", "CRISPR Guide Capture",
    "Antigen Capture",
})


def run_mux_tests(var_names, class_df, stats_df, obsm_col_sets, ref_map,
                  ref_provenance, full_distro=True, bookkeeping=None,
                  var_ids_list=None, unsupervised=None):
    """T-MUX-01 (multiplexing-tag placement), T-MUX-02 (demux bookkeeping),
    T-MUX-03 (unsupervised counts-only discovery), T-REF-09 (reconciliation).

    A multiplexing tag is identified by UNION of independent arms (no single
    arm is trusted alone for custom names):
      R  ref/var type == 'Multiplexing Capture' (any name — custom-safe);
      ID feature-ID column matching HTO_ID_RE (TotalSeq HTO IDs, Hashtag IDs;
         catches HTOs labeled 'Antibody Capture', e.g. TotalSeq-B);
      S  hto-like spelling of the var_name (weak fallback);
      B  bookkeeping attestation: uns hto-list values or obs classification
         labels that occur in var/obsm.
    obsm_col_sets: {obsm_key: [columns]} for obsm slots holding DataFrames.
    stats_df must cover at least all adt-like/hto-like var_names plus all
    non-GEX reference names in var (minimal fallback) or all features (full
    distro). bookkeeping: summarize_mux_bookkeeping() output or None.
    """
    tests = []
    stat = {r.var_name: r for r in stats_df.itertuples()}
    obsm_all = {c for cols in obsm_col_sets.values() for c in cols}
    var_set = set(var_names)

    ref_mux = sorted({n for n, t in (ref_map or {}).items()
                      if t == "Multiplexing Capture"})
    spell_mux = sorted({n for n in var_names if "hto-like" in tag_name(n)[1]
                        and n not in (ref_map or {})})
    id_mux = sorted({n for n, i in zip(var_names, var_ids_list or [])
                     if i and HTO_ID_RE.search(str(i))})
    bk = bookkeeping or {}
    bk_tags = sorted(({t for v in (bk.get("uns_hto_lists") or {}).values()
                       for t in v}
                      | {v for v in (bk.get("class_values") or [])})
                     & (var_set | obsm_all))
    mux_names = sorted(set(ref_mux) | set(spell_mux) | set(id_mux) | set(bk_tags))

    def _counts(names):
        rows = [(n, stat[n].mean) for n in names if n in stat]
        return ([n for n, m in rows if m and m > 0],
                [n for n, m in rows if not m])

    mux_nonzero, mux_zero = _counts(mux_names)
    mux_in_obsm = sorted({c for cols in obsm_col_sets.values() for c in cols}
                         & set(mux_names))
    src_note = (f"ref:{len(ref_mux)} htoID:{len(id_mux)} spelling:{len(spell_mux)} "
                f"bookkeeping:{len(bk_tags)}"
                + ("" if ref_map else " (NO feature reference: custom "
                   "patient/sample-ID hashtags are undiscoverable — "
                   "spelling scan only)"))

    # T-MUX-01: multiplexing-tag placement
    if mux_nonzero:
        tests.append(("T-MUX-01", "MUX", "FAIL",
                      f"{len(mux_nonzero)}/{len(mux_names)} multiplexing-tag names "
                      f"({src_note}) have NONZERO counts in X — possible HTO/GEX "
                      "mixing; curator review required, nothing removed",
                      len(mux_nonzero), mux_nonzero[:10]))
    elif mux_in_obsm:
        tests.append(("T-MUX-01", "MUX", "PASS",
                      f"{len(mux_names)} multiplexing-tag names ({src_note}) "
                      f"resolve to obsm ({mux_in_obsm[:5]}); "
                      "none with counts in X",
                      0, mux_in_obsm[:10]))
    elif mux_zero:
        tests.append(("T-MUX-01", "MUX", "PASS",
                      f"{len(mux_zero)} multiplexing-tag names ({src_note}) present "
                      "in X/var but ALL zero-count (empty placeholder columns); "
                      "no mixing evidence; retained",
                      0, mux_zero[:10]))
    elif ref_map is None:
        tests.append(("T-MUX-01", "MUX", "PASS",
                      "no multiplexing tags discovered and no feature reference "
                      "file available: custom-named hashtags cannot be ruled in "
                      "or out from spelling alone; HTO spelling scan clean",
                      0, []))
    else:
        tests.append(("T-MUX-01", "MUX", "PASS",
                      f"feature reference lists no Multiplexing Capture rows "
                      f"({ref_provenance}); HTO spelling scan clean",
                      0, []))

    # T-REF-09: feature-barcode reconciliation — every non-GEX reference row
    # must be accounted for (obsm, zero-count var, or nonzero var).
    if ref_map:
        ref_bc = sorted({n for n, t in ref_map.items() if t in NON_GEX_FEATURE_TYPES})
        obsm_all = {c for cols in obsm_col_sets.values() for c in cols}
        var_set = set(var_names)
        in_obsm = [n for n in ref_bc if n in obsm_all]
        in_var_nz = [n for n in ref_bc if n in var_set
                     and n in stat and stat[n].mean and stat[n].mean > 0]
        in_var_zero = [n for n in ref_bc if n in var_set and n not in in_var_nz]
        missing = [n for n in ref_bc if n not in obsm_all and n not in var_set]
        accounted = len(in_obsm) + len(in_var_nz) + len(in_var_zero)
        tests.append(("T-REF-09", "REF",
                      "WARN" if missing else "PASS",
                      f"feature-barcode reconciliation ({ref_provenance}): "
                      f"{len(ref_bc)} non-GEX reference rows; {len(in_obsm)} in obsm, "
                      f"{len(in_var_nz)} nonzero in X/var, {len(in_var_zero)} zero-count "
                      f"in X/var, {len(missing)} MISSING everywhere (possible data loss)"
                      if missing else
                      f"feature-barcode reconciliation ({ref_provenance}): all "
                      f"{len(ref_bc)} non-GEX reference rows accounted for "
                      f"({len(in_obsm)} obsm, {len(in_var_nz)} nonzero X/var, "
                      f"{len(in_var_zero)} zero X/var)",
                      len(missing), missing[:10]))
    else:
        tests.append(("T-REF-09", "REF", "PASS",
                      "not applicable (no feature reference discovered)", 0, []))

    # T-REF-10: gene-name / feature-type conflicts. A CITE/HTO tag carrying a
    # genuine gene symbol or ENSEMBL ID classifies as a confirmed transcript,
    # so without this test the exemption would hide it and the counts would
    # silently corrupt that gene's apparent expression. Two checks:
    # (a) confirmed transcript AND declared non-GEX barcode: nonzero X counts
    #     -> FAIL (silent corruption); zero counts -> WARN (placeholder);
    #     counts unavailable -> WARN (unchecked).
    # (b) exact var_names/obsm-column collisions across modalities (e.g. a tag
    #     counted in both X and an antibody matrix) -> WARN (ambiguous identity).
    conf = class_df[class_df["confirmed"]
                    & class_df["ref_feature_type"].isin(NON_GEX_FEATURE_TYPES)]
    conf_nz, conf_zero, conf_unchecked = [], [], []
    for r in conf.itertuples():
        if r.var_name in stat and stat[r.var_name].mean \
                and stat[r.var_name].mean > 0:
            conf_nz.append(f"{r.var_name}[{r.ref_feature_type},"
                           f"mean={stat[r.var_name].mean:.1f}]")
        elif r.var_name in stat:
            conf_zero.append(r.var_name)
        else:
            conf_unchecked.append(r.var_name)
    overlap = sorted(var_set & obsm_all)
    problems = conf_nz + conf_zero + conf_unchecked + overlap
    if conf_nz:
        tests.append(("T-REF-10", "REF", "FAIL",
                      f"{len(conf_nz)} confirmed gene transcripts are declared "
                      "non-GEX feature barcodes yet carry NONZERO counts in X — "
                      "silent corruption of those genes' expression; curator "
                      "review required, nothing removed",
                      len(problems), (conf_nz + overlap)[:10]))
    elif problems:
        tests.append(("T-REF-10", "REF", "WARN",
                      "gene-name/feature-type conflicts "
                      f"(zero-count:{conf_zero[:5]} unchecked:{conf_unchecked[:5]}; "
                      f"var/obsm collisions:{overlap[:5]}); retained",
                      len(problems), problems[:10]))
    else:
        tests.append(("T-REF-10", "REF", "PASS",
                      "no gene-name/feature-type conflicts; "
                      "no var/obsm name collisions",
                      0, []))

    # T-MUX-02: demux bookkeeping consistency — the traces a hashtag workflow
    # leaves in obs/uns must agree with the discovered tags. Checks (WARN =
    # curator review, never a filter): classification labels unknown anywhere;
    # uns-listed tags with counts nowhere; filtered-listed tags still broadly
    # present (frac > 0.5 review trigger).
    bk = bookkeeping or {}
    bk_class_vals = bk.get("class_values") or []
    bk_uns_hto = bk.get("uns_hto_lists") or {}
    bk_uns_filt = bk.get("uns_filtered_lists") or {}
    bk_overlap = bk.get("tag_obs_var_overlap") or []
    if not bk or (not bk_class_vals and not bk_uns_hto and not bk_uns_filt
                  and not bk_overlap):
        tests.append(("T-MUX-02", "MUX", "PASS",
                      "no demux bookkeeping detected in obs/uns "
                      "(no classification column, tag lists, or obs/var tag overlap)",
                      0, []))
    else:
        universe = (set(ref_mux) | set(spell_mux) | set(mux_in_obsm)
                    | {t for v in bk_uns_hto.values() for t in v}
                    | {t for v in bk_uns_filt.values() for t in v}
                    | set(bk_overlap))
        obsm_all = {c for cols in obsm_col_sets.values() for c in cols}
        var_set = set(var_names)
        unknown_labels = sorted({v for v in bk_class_vals
                                 if v not in universe
                                 and v.lower() not in ("doublet", "negative", "unassigned",
                                                       "ambiguous", "multiplet")})
        listed_nowhere = sorted({t for v in bk_uns_hto.values() for t in v}
                                - obsm_all - var_set)
        filt_broad, filt_nocount = [], []
        for v in bk_uns_filt.values():
            for t in v:
                if t in stat and np.isfinite(stat[t].frac_expressed):
                    if stat[t].frac_expressed > 0.5:
                        filt_broad.append(f"{t}(frac={stat[t].frac_expressed:.2f})")
                elif t in var_set or t in obsm_all:
                    filt_nocount.append(f"{t}(counts unchecked in minimal mode)")
        problems = (unknown_labels + listed_nowhere + filt_broad + filt_nocount)
        detail = (f"class_col={bk.get('class_column')}:{len(bk_class_vals)} labels; "
                  f"uns kept lists={ {k: len(v) for k, v in bk_uns_hto.items()} }; "
                  f"uns filtered lists={ {k: len(v) for k, v in bk_uns_filt.items()} }; "
                  f"obs/var tag overlap={len(bk_overlap)}")
        if problems:
            tests.append(("T-MUX-02", "MUX", "WARN",
                          f"demux bookkeeping inconsistencies ({detail}): "
                          f"unknown labels={unknown_labels[:5]}; "
                          f"listed-but-countless={listed_nowhere[:5]}; "
                          f"filtered-but-broad={filt_broad[:5]}"
                          + ("; counts unchecked for "
                             f"{filt_nocount[:5]}" if filt_nocount else "")
                          + " (retained)",
                          len(problems), problems[:10]))
        else:
            tests.append(("T-MUX-02", "MUX", "PASS",
                          f"demux bookkeeping consistent ({detail})",
                          0, []))

    # T-MUX-03: unsupervised counts-only discovery. WARN-only by design:
    # statistical suspicion is never identity evidence, so this test can only
    # point, never convict (and never FAIL).
    if unsupervised is None:
        tests.append(("T-MUX-03", "MUX", "PASS",
                      "unsupervised discovery skipped (--no-distro)", 0, []))
        uns_panel_note = ""
    elif unsupervised.get("decision") == "panel":
        ev_of = dict(zip(class_df["var_name"], class_df["evidence"]))
        bio_of = dict(zip(class_df["var_name"], class_df["biotype"]))
        conf_of = dict(zip(class_df["var_name"], class_df["confirmed"]))
        panel = unsupervised["panel"]
        overlap = sorted(set(panel) & set(mux_names))
        def _bio_disp(nm):
            b = bio_of.get(nm)
            return "no-biotype" if _missing(b) else b

        ex = [f"{c['var_name']}[group={c['group_n']},"
              f"{ev_display(ev_of.get(c['var_name'], '?'))},"
              f"{_bio_disp(c['var_name'])}"
              f"{',CONFIRMED-GENE!' if conf_of.get(c['var_name']) else ''}]"
              for c in unsupervised["candidates"][:8]]
        tests.append(("T-MUX-03", "MUX", "WARN",
                      f"unsupervised counts-only discovery points at "
                      f"{len(panel)} hashtag-like columns "
                      f"(partition: {unsupervised['reason']}; "
                      f"{len(overlap)}/{len(panel)} overlap the authoritative "
                      "mux set — independent confirmation, not proof). "
                      "Review list, retained; identity requires reference or "
                      "bookkeeping evidence",
                      len(panel), ex))
        uns_panel_note = (f" unsupervised discovery independently points at "
                          f"{len(panel)}: {', '.join(panel[:8])}.")
    else:
        tests.append(("T-MUX-03", "MUX", "PASS",
                      "unsupervised counts-only discovery: "
                      f"{unsupervised.get('reason', 'no panel')}",
                      0, []))
        uns_panel_note = ""

    # MUX placement declaration
    if mux_names:
        placement = (f"multiplexing tags discovered ({len(ref_mux)} ref "
                     f"Multiplexing Capture, {len(id_mux)} HTO-ID, "
                     f"{len(spell_mux)} spelling, {len(bk_tags)} bookkeeping): "
                     f"{len(mux_nonzero)} nonzero in X/var, "
                     f"{len(mux_zero)} zero-count in X/var, "
                     f"{len(mux_in_obsm)} in obsm." + uns_panel_note)
    elif ref_map is not None:
        placement = ("no multiplexing tags: feature reference lists no "
                     "Multiplexing Capture rows and HTO spelling scan is clean."
                     + uns_panel_note)
    else:
        placement = ("no multiplexing tags discovered and no feature reference "
                     "available: custom-named hashtags (patient/sample IDs) "
                     "cannot be excluded by spelling alone." + uns_panel_note)
    return tests, placement


# ============================================================================
# Unsupervised hashtag discovery (counts only — no names, references, or bookkeeping)
# Hashtag panels have a distinctive joint signature: each tag is bimodal
# (ambient background + bright positive mode) AND the tags are mutually
# exclusive per cell (winner-takes-all: each cell's argmax is its hashtag,
# with a large top1/top2 margin). Genes do not do this jointly — a globally
# high gene wins argmax nearly everywhere (degenerate partition) with ~1x
# margins. Output is a RANKED REVIEW LIST (WARN-only downstream): statistical
# suspicion is never identity evidence, so discovery alone can never FAIL.
# ============================================================================

BIMODAL_BC_CUTOFF = 0.555  # standard bimodality-coefficient threshold


def discover_hashtag_candidates(stats_df, X, var_names, max_preselect=200,
                                min_group_frac=0.01, margin_p10_min=5.0,
                                winner_median_min=20.0):
    """Point at hashtag-like columns using counts alone.

    Inputs used: per-feature stats (max/mean/bc) + the count matrix Itself.
    Deliberately NOT used: var_names spellings, var type/ID columns, obs, uns,
    references. Returns a dict:
      searched (bool), n_preselect, cutoff_max, candidates (list of dicts with
      var_name/group_n/margin stats per candidate... group-level),
      panel (ranked names), margin_median, margin_p10, decision, reason.
    decision is 'panel' or 'none', with a human-readable reason.
    """
    out = {"searched": True, "n_preselect": 0, "cutoff_max": None,
           "candidates": [], "panel": [], "margin_median": None,
           "margin_p10": None, "decision": "none", "reason": ""}
    nz = stats_df[stats_df["mean"] > 0]
    if not len(nz):
        out["reason"] = "no nonzero features"
        return out
    med_max = float(nz["max"].median())
    cutoff = max(50.0, 10.0 * med_max)
    out["cutoff_max"] = cutoff
    pre = nz[(nz["max"] >= cutoff) & np.isfinite(nz["bc"])
             & (nz["bc"] > BIMODAL_BC_CUTOFF)]
    pre = pre.nlargest(max_preselect, "max")
    out["n_preselect"] = int(len(pre))
    if len(pre) < 2:
        out["reason"] = (f"only {len(pre)} bimodal high-max columns "
                         f"(cutoff max>={cutoff:.0f}, BC>{BIMODAL_BC_CUTOFF})")
        return out

    # Winner-takes-all partition test on preselected columns only.
    # Columns may be scattered: gather via contiguous runs (never the span).
    name_to_j = {n: j for j, n in enumerate(var_names)}
    idx = [name_to_j[n] for n in pre["var_name"] if n in name_to_j]
    if len(idx) < 2:
        out["reason"] = "preselected names not aligned to X"
        return out
    ordered = sorted(set(idx))
    blocks = []
    rs = re_ = ordered[0]
    for j in ordered[1:]:
        if j == re_ + 1:
            re_ = j
        else:
            blocks.append((rs, re_))
            rs = re_ = j
    blocks.append((rs, re_))
    offsets, pos = {}, 0
    parts = []
    for a, b in blocks:
        parts.append(_dense_slice(X, a, b + 1))
        for j in range(a, b + 1):
            offsets[j] = pos
            pos += 1
    wide = np.concatenate(parts, axis=1)
    M = np.column_stack([wide[:, offsets[j]] for j in idx])
    n_cells = M.shape[0]
    min_group = max(10, int(min_group_frac * n_cells))

    alive = list(range(M.shape[1]))
    while len(alive) >= 2:
        top = M[:, alive].argmax(axis=1)
        counts = np.bincount(top, minlength=len(alive))
        keep = [alive[k] for k in range(len(alive)) if counts[k] >= min_group]
        if len(keep) == len(alive):
            break
        alive = keep
    if len(alive) < 2:
        out["reason"] = (f"partition degenerate: fewer than 2 columns win "
                         f">={min_group} cells each")
        return out

    top = M[:, alive].argmax(axis=1)
    srt = np.sort(M[:, alive], axis=1)
    margin = (srt[:, -1] + 1.0) / (srt[:, -2] + 1.0)
    med_margin = float(np.median(margin))
    p10_margin = float(np.percentile(margin, 10))
    out["margin_median"] = med_margin
    out["margin_p10"] = p10_margin
    counts = np.bincount(top, minlength=len(alive))
    ranked = sorted(zip([pre["var_name"].iloc[k] for k in alive], counts),
                    key=lambda t: -t[1])
    # Winner-count floor (feature-barcode scale prior): oligo-derived tags
    # yield order-of-magnitude higher UMI counts than transcripts, so the
    # median winning count in every group must clear it. This separates
    # hashtag panels from sex-split genes (XIST/RPS4Y1 win groups at
    # single-digit counts).
    # NOTE: `top` holds positions into `alive` (argmax over M[:, alive]).
    pos_of = {}
    for k, kk in enumerate(alive):
        pos_of[pre["var_name"].iloc[kk]] = k
    winner_med = {}
    for n, _ in ranked:
        k = pos_of[n]
        cells = top == k
        winner_med[n] = float(np.median(M[cells][:, alive[k]]))
    out["winner_medians"] = {n: round(m, 1) for n, m in winner_med.items()}
    min_winner = min(winner_med.values())
    out["candidates"] = [{"var_name": n, "group_n": int(c),
                          "group_frac": round(float(c) / n_cells, 4),
                          "winner_median": out["winner_medians"][n]}
                         for n, c in ranked]
    out["panel"] = [n for n, _ in ranked]
    if p10_margin < margin_p10_min:
        out["reason"] = (f"panel of {len(ranked)} splits cells but the winner's "
                         f"lead is weak (p10 top1/top2={p10_margin:.1f}x < "
                         f"{margin_p10_min}x): not hashtag-like")
        return out
    if min_winner < winner_median_min:
        out["reason"] = (f"panel of {len(ranked)} splits cells with a strong lead "
                         f"but winners are transcript-scale (min group winner "
                         f"median={min_winner:.1f} < {winner_median_min}): "
                         "not hashtag-like")
        return out
    out["decision"] = "panel"
    out["reason"] = (f"panel of {len(ranked)} columns partitions cells into "
                     f"sizable disjoint groups (min group "
                     f"{min(int(c) for _, c in ranked)} cells) with a strong "
                     f"winner lead (median top1/top2={med_margin:.0f}x, "
                     f"p10={p10_margin:.0f}x)")
    return out


# ============================================================================
# Distribution statistics (read-only; backed-safe)
# ============================================================================

def _dense_slice(X, lo, hi):
    """Return X[:, lo:hi] as a dense float64 ndarray.

    Handles backed h5ad _CSRDataset (which scipy.sparse.issparse does not
    recognise) by calling .toarray() directly.
    """
    sl = X[:, lo:hi]
    if hasattr(sl, "toarray"):
        return np.asarray(sl.toarray(), dtype=np.float64)
    arr = np.asarray(sl, dtype=np.float64)
    if arr.ndim == 0:  # 0-dim object array guard — must not silently pass
        raise ValueError(
            f"Dense slice [{lo}:{hi}] produced a 0-dim array; refusing to "
            "continue rather than report wrong statistics."
        )
    return arr


def compute_distribution_stats(X, var_names, chunk_size=500):
    """Per-feature mean/var/CV/overdispersion/fraction-expressed/max/BC.

    BC is the bimodality coefficient (skew^2+1)/kurtosis (Pearson); values
    above ~0.555 suggest bimodal/multimodal counts. NaN where undefined
    (e.g. zero-count columns). Raises on unexpected shapes instead of guessing.
    """
    from scipy.stats import kurtosis, skew
    n_cells, n_genes = X.shape
    if n_genes != len(var_names):
        raise ValueError(
            f"X has {n_genes} columns but {len(var_names)} var_names; "
            "refusing to compute misaligned statistics."
        )
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    # Keep peak chunk memory bounded for large cell counts.
    if n_cells > 100_000 and chunk_size > 200:
        chunk_size = 200

    records = []
    for lo in range(0, n_genes, chunk_size):
        hi = min(lo + chunk_size, n_genes)
        chunk = _dense_slice(X, lo, hi)
        if chunk.shape != (n_cells, hi - lo):
            raise ValueError(
                f"Slice [{lo}:{hi}] has shape {chunk.shape}, expected "
                f"({n_cells}, {hi - lo}); refusing to continue."
            )
        means = chunk.mean(axis=0)
        variances = chunk.var(axis=0)
        frac = (chunk > 0).mean(axis=0)
        maxs = chunk.max(axis=0)
        with np.errstate(all="ignore"):
            sk = skew(chunk, axis=0)
            ku = kurtosis(chunk, axis=0, fisher=False)  # Pearson kurtosis
            bc = (sk ** 2 + 1) / np.where(ku != 0, ku, np.nan)
        for j in range(chunk.shape[1]):
            mean = float(means[j])
            var = float(variances[j])
            records.append({
                "var_name": var_names[lo + j],
                "mean": mean,
                "var": var,
                "cv": float(np.sqrt(var) / mean) if mean > 0 else 0.0,
                "overdisp": float(var / mean) if mean > 0 else 0.0,
                "frac_expressed": float(frac[j]),
                "max": float(maxs[j]),
                "bc": float(bc[j]) if np.isfinite(bc[j]) else float("nan"),
            })
    return pd.DataFrame(records)


# ============================================================================
# Test suite: REF (references) and DIST (distributions)
# Each test returns (result, summary, n_affected, examples).
# Results: PASS / WARN / FAIL. FAIL = structural corruption or ADT mixing
# evidence. WARN = needs curator review (incl. expected reference-version drift).
# ============================================================================

def run_ref_tests(var_names, class_df, var_ensg_list, ensg_col, gencode_ensg,
                  gmaps):
    tests = []

    # T-REF-01: uniqueness
    dup_mask = class_df["var_name"].duplicated(keep=False)
    n_dup = int(dup_mask.sum())
    tests.append(("T-REF-01", "REF",
                  "FAIL" if n_dup else "PASS",
                  f"{n_dup} duplicated var_names (duplicates corrupt X/var alignment)"
                  if n_dup else "var_names are unique",
                  n_dup,
                  sorted(class_df[dup_mask]["var_name"].unique().tolist()[:10])))

    # T-REF-02: empty / blank names
    blanks = [n for n in var_names if not str(n).strip()]
    tests.append(("T-REF-02", "REF",
                  "FAIL" if blanks else "PASS",
                  f"{len(blanks)} empty/blank var_names" if blanks
                  else "no empty/blank var_names",
                  len(blanks), blanks[:10]))

    # T-REF-03: malformed ENSG-like names
    malformed = [n for n in var_names
                 if ENSG_LIKE_RE.match(n) and not ENSG_STRICT_RE.match(n)]
    tests.append(("T-REF-03", "REF",
                  "FAIL" if malformed else "PASS",
                  f"{len(malformed)} ENSG-like names with invalid format"
                  if malformed else "all ENSG-like names well-formed",
                  len(malformed), malformed[:10]))

    # T-REF-04: format-valid ENSG IDs absent from GENCODE v44 (deprecated/newer)
    all_ensg = [e for e in (var_ensg_list or []) if e]
    unknown_ensg = sorted({e for e in all_ensg if e not in gencode_ensg})
    tests.append(("T-REF-04", "REF",
                  "WARN" if unknown_ensg else "PASS",
                  f"{len(unknown_ensg)} distinct ENSG IDs not in GENCODE v44 "
                  "(deprecated or newer than v44; retained)"
                  if unknown_ensg else "all observed ENSG IDs are in GENCODE v44",
                  len(unknown_ensg), unknown_ensg[:10]))

    # T-REF-05: symbol<->ENSG cross-column consistency (only where var has both)
    if ensg_col and var_ensg_list:
        e2s = gmaps.get("ensg_to_symbol", {})
        tested, no_proper, renamed, mism = 0, [], [], []
        for row in class_df.itertuples():
            if row.ensg and row.ensg in e2s:
                tested += 1
                v44sym = e2s[row.ensg]
                if row.var_name != v44sym:
                    m = f"{row.var_name}!=GENCODE:{v44sym}({row.ensg})"
                    mism.append(m)
                    (no_proper if v44sym == row.ensg else renamed).append(m)
        tests.append(("T-REF-05", "REF",
                      "WARN" if mism else "PASS",
                      f"{len(mism)}/{tested} symbol/ENSG pairs disagree with "
                      "GENCODE v44 (retained): "
                      f"{len(no_proper)} where v44 itself has no proper symbol "
                      "(ENSG-as-name; study's older symbol is more informative), "
                      f"{len(renamed)} where v44 renamed the locus. "
                      "Expected under reference-version drift, e.g. GRCh37-era "
                      "annotations"
                      if mism else
                      (f"all {tested} symbol/ENSG pairs agree with GENCODE v44"
                       if tested else "no ENSG-anchored pairs to check"),
                      len(mism), mism[:10]))
    else:
        tests.append(("T-REF-05", "REF", "PASS",
                      "not applicable (no ENSG var column)", 0, []))

    # T-REF-06: ncRNA authority check — ncRNA-pattern names vs GENCODE biotypes
    ncrna_rows = class_df[class_df["primary_tag"].isin(NCRNA_PATTERN_TAGS)]
    confirmed, other_bt, unconfirmed = [], [], []
    for row in ncrna_rows.itertuples():
        bt_missing = _missing(row.biotype)
        bts = set() if bt_missing else set(str(row.biotype).split("|"))
        if bts & set(NCRNA_BIOTYPES):
            confirmed.append(row.var_name)
        elif not bt_missing:
            other_bt.append(f"{row.var_name}[{row.biotype}]")
        else:
            unconfirmed.append(row.var_name)
    tests.append(("T-REF-06", "REF",
                  "WARN" if unconfirmed else "PASS",
                  f"ncRNA-pattern names: {len(confirmed)} confirmed ncRNA biotype, "
                  f"{len(other_bt)} other GENCODE biotype, {len(unconfirmed)} "
                  "unconfirmed (absent from GENCODE v44; retained)"
                  if len(ncrna_rows) else "no ncRNA-pattern names present",
                  len(unconfirmed),
                  (unconfirmed[:10] + [f"other-biotype:{x}" for x in other_bt[:5]])))

    # T-REF-07: GENCODE "artifact" biotype — flagged by GENCODE itself (19 genes)
    art = class_df[class_df["biotype"] == "artifact"]
    tests.append(("T-REF-07", "REF",
                  "WARN" if len(art) else "PASS",
                  f"{len(art)} features resolve to GENCODE gene_type 'artifact' "
                  "(curator review; retained)"
                  if len(art) else "no features with GENCODE 'artifact' biotype",
                  int(len(art)), art["var_name"].head(10).tolist()))

    # T-REF-08: authority coverage — guards against corrupt builds
    # (e.g. wrong column used as var_names would strand nearly all in E5)
    resolved = int(class_df["evidence"].isin(
        ["E1_ensg_match", "E2_gencode_symbol", "E3a_hgnc_renamed",
         "E3b_hgnc_current", "E4_ncbi_alias"]).sum())
    rate = resolved / max(len(class_df), 1)
    tests.append(("T-REF-08", "REF",
                  "WARN" if rate < 0.50 else "PASS",
                  f"authority-resolved {resolved}/{len(class_df)} "
                  f"({100 * rate:.1f}%); unresolved remainder retained with tags",
                  len(class_df) - resolved, []))

    return tests


def run_dist_tests(class_df, stats_df, adt_like_names, adt_panel,
                   full_distro=True, confirmed_names=None, panel_names=None):
    """T-DIST-01 (envelope), T-DIST-02/03/04.

    confirmed_names: var_names with authoritative gene identity (the
    `confirmed` classification flag). They are exempt from suspicion as
    HTO/ADT *suspects* — except corroborated distribution evidence: a
    confirmed transcript that BOTH falls inside the protein envelope AND
    joins the unsupervised hashtag panel carries a distribution that
    contradicts its gene identity, and the curator must be made aware
    (WARN, never concealed). Uncorroborated housekeeping overlap is reported
    as a count only. Globals (T-DIST-03/04) always use all features.
    panel_names: unsupervised hashtag panel (set) or empty.
    """
    tests = []
    stat = {r.var_name: r for r in stats_df.itertuples()}

    def feats(names):
        out = []
        for n in names:
            s = stat.get(n)
            if s is not None:
                out.append((n, s.mean, s.overdisp, s.frac_expressed, s.bc,
                            s.max))
        return out

    # T-DIST-01: envelope overlap vs the REAL ADT panel (calibrated, not invented)
    if not full_distro:
        tests.append(("T-DIST-01", "DIST", "PASS",
                      "skipped (--no-distro); rerun without the flag for the "
                      "empirical-envelope test", 0, []))
        inside = []
    elif adt_panel is None:
        tests.append(("T-DIST-01", "DIST", "WARN",
                      "ADT reference panel missing (cache/adt_reference_panel.json); "
                      "envelope test skipped", 0, []))
        inside = []
    else:
        env = adt_panel["envelope_p5_p95"]
        confirmed = confirmed_names or set()
        panel = panel_names or set()

        def within(v, lo, hi):
            return np.isfinite(v) and lo <= v <= hi

        def _inside(n):
            s = stat[n]
            return (within(s.mean, env["mean"]["p5"], env["mean"]["p95"])
                    and within(s.overdisp, env["overdisp"]["p5"], env["overdisp"]["p95"])
                    and within(s.frac_expressed, env["frac_expressed"]["p5"],
                               env["frac_expressed"]["p95"])
                    and within(s.bc, env["bc"]["p5"], env["bc"]["p95"]))

        # Suspects: non-confirmed features inside the envelope.
        testable = [n for n in stats_df["var_name"] if n not in confirmed]
        n_exempt = len(stats_df) - len(testable)
        inside = [n for n in testable if _inside(n)]
        ev_of = dict(zip(class_df["var_name"], class_df["evidence"]))
        tier_counts = pd.Series([ev_display(ev_of.get(n, "?"))
                                 for n in inside]).value_counts().to_dict()
        e5_inside_nonzero = [n for n in inside
                             if ev_of.get(n) == "E5_unresolved_retained"
                             and stat[n].mean > 0]
        # Awareness: confirmed transcripts are never suspects, BUT a confirmed
        # gene that is both protein-distributed AND partitions with a hashtag
        # panel contradicts its own identity — the curator must see it.
        confirmed_inside = sorted(n for n in stats_df["var_name"]
                                  if n in confirmed and _inside(n))
        corroborated = sorted(set(confirmed_inside) & set(panel))
        plain_overlap = len(confirmed_inside) - len(corroborated)
        parts = [f"{len(inside)}/{len(testable)} non-confirmed features fall "
                 f"inside the empirical ADT envelope {adt_panel['provenance']}; "
                 f"{n_exempt} confirmed transcripts exempted (known gene "
                 f"identity, never HTO/ADT suspects); tiers={tier_counts}"]
        if e5_inside_nonzero:
            parts.append(f"{len(e5_inside_nonzero)} unresolved inside with "
                         "nonzero counts -> curator review list (retained)")
        if corroborated:
            parts.append(f"{len(corroborated)} CONFIRMED transcripts inside the "
                         "envelope AND in the unsupervised hashtag panel "
                         f"({', '.join(corroborated[:8])}): protein-like "
                         "distribution on a gene identity — possible gene-named "
                         "tag, curator must rule it out (retained)")
        elif confirmed_inside:
            parts.append(f"+{plain_overlap} confirmed inside with no panel "
                         "corroboration (expected expression overlap)")
        tests.append(("T-DIST-01", "DIST",
                      "WARN" if (e5_inside_nonzero or corroborated) else "PASS",
                      "; ".join(parts),
                      len(e5_inside_nonzero) + len(corroborated),
                      (e5_inside_nonzero + corroborated)[:10]))

    # T-DIST-02: ADT-like spelling AND nonzero counts = mixing evidence
    rows = feats(adt_like_names)
    nonzero = [n for n, mean, *_ in rows if mean and mean > 0]
    zero = [n for n, mean, *_ in rows if not mean]
    if not adt_like_names:
        tests.append(("T-DIST-02", "DIST", "PASS",
                      "no ADT-like var_names in X", 0, []))
    elif nonzero:
        tests.append(("T-DIST-02", "DIST", "FAIL",
                      f"{len(nonzero)}/{len(rows)} ADT-like var_names have NONZERO "
                      "counts in X — possible ADT/GEX mixing; curator review required, "
                      "nothing removed by this script",
                      len(nonzero), nonzero[:10]))
    else:
        tests.append(("T-DIST-02", "DIST", "PASS",
                      f"{len(zero)} ADT-like var_names present but ALL zero-count "
                      "(empty placeholder columns, consistent with feature-reference "
                      "carryover); no mixing evidence; retained",
                      0, zero[:10]))

    # T-DIST-03: global sanity — GEX sparsity vs ADT sparsity
    if not full_distro:
        tests.append(("T-DIST-03", "DIST", "PASS",
                      "skipped (--no-distro)", 0, []))
    elif adt_panel is not None:
        panel_fracs = [f["frac_expressed"] for f in adt_panel["features"]]
        adt_med_frac = float(np.median(panel_fracs))
        gex_med_frac = float(stats_df["frac_expressed"].median())
        gex_med_mean = float(stats_df["mean"].median())
        tests.append(("T-DIST-03", "DIST", "PASS",
                      f"median frac_expressed: GEX={gex_med_frac:.3f} vs ADT panel={adt_med_frac:.3f}; "
                      f"median GEX mean={gex_med_mean:.3f} (ADT panel median mean="
                      f"{np.median([f['mean'] for f in adt_panel['features']]):.2f}). "
                      "GEX is far sparser, as expected for RNA vs protein counts",
                      0, []))
    else:
        tests.append(("T-DIST-03", "DIST", "WARN",
                      "ADT panel missing; global comparison skipped", 0, []))

    # T-DIST-04: zero-count columns (informational; QC handles them)
    n_zero = int((stats_df["mean"] == 0).sum())
    tests.append(("T-DIST-04", "DIST", "PASS",
                  f"{n_zero}/{len(stats_df)} columns are all-zero "
                  "(no signal; naturally excluded by any variance filter; retained)",
                  n_zero, stats_df[stats_df["mean"] == 0]["var_name"].head(5).tolist()
                  if n_zero else []))

    # ADT placement verdict string (storage locations)
    if not adt_like_names:
        verdict = ("No ADT-like var_names in the GEX matrix. "
                   "No evidence of ADT counts mixed into X.")
    elif nonzero:
        verdict = (f"{len(nonzero)}/{len(rows)} ADT-like var_names have NONZERO counts "
                   "in X — FAIL (antibody counts in RNA matrix). Human review required; nothing removed.")
    else:
        verdict = (f"{len(zero)} ADT-like var_names present but ALL zero counts "
                   "(empty placeholders). No evidence of ADT counts mixed into GEX. "
                   "Retained as-is.")
    return tests, verdict


# ============================================================================
# Per-dataset driver
# ============================================================================

def verify_dataset(ds, base_dir, gencode_symbols, gencode_ensg, gmaps,
                   hgnc_cache, ncbi_cache, adt_panel,
                   chunk_size=500, run_distro=True,
                   ref_map=None, ref_provenance="", h5ad_path=None,
                   verbose=False, brief=False):
    """Run the full verification on one dataset h5ad (read-only, backed).

    ref_map ({feature_name: feature_type}) comes from discover_feature_reference
    and is merged UNDER the in-matrix var feature_types column (the var column
    describes exactly what is in THIS X, including custom hashtag names).
    h5ad_path optionally points directly at a file (label defaults to ds or
    the file stem); otherwise file_database/<DS>/<DS>.h5ad is used.
    verbose=True prints full tables; default prints a compact readable digest
    (full tables always land in the CSVs via --out).
    Returns a dict with n_cells, n_genes, classification, standard_table,
    test_table, adt_placement, mux_placement, adt_verdict and overall.
    Returns None if the h5ad is missing.
    """
    import builtins
    _print = builtins.print
    if brief:
        # Silence the guided tour below (it still computes everything the
        # verdicts need); header, mini-inventory and verdicts use _print.
        print = lambda *args, **kwargs: None  # noqa: silence tour output
    else:
        print = _print
    if h5ad_path is None:
        h5ad_path = os.path.join(base_dir, ds, f"{ds}.h5ad")
        label = ds
    else:
        label = ds or os.path.splitext(os.path.basename(h5ad_path))[0]
    _print(f"\n{'=' * 70}")
    _print(f"  {label}")
    _print(f"{'=' * 70}")

    if not os.path.exists(h5ad_path):
        _print("  SKIP: h5ad not found (no data touched)")
        return None

    import anndata as ad
    adata = ad.read_h5ad(h5ad_path, backed="r")  # read-only
    try:
        var_names = adata.var_names.tolist()
        n_cells, n_genes = adata.X.shape
        _print(f"  {n_cells} cells x {n_genes} genes (read-only, backed)")

        if brief:
            obsm_keys = list(adata.obsm.keys()) if hasattr(adata, "obsm") else []
            adt_obsm = [k for k in obsm_keys
                        if any(s in k.lower() for s in ("adt", "antibody", "protein"))]
            if "feature_types" in list(adata.var.columns):
                _print(f"  var feature_types: "
                       f"{adata.var['feature_types'].value_counts().to_dict()}")
        else:
            # --- Slot inventory: what attribute types does this dataset carry? ---
            print("\n  [Slots] AnnData attribute inventory")

            def _short(items, k=8):
                items = list(items)
                head = ", ".join(map(str, items[:k]))
                tail = f" +{len(items) - k} more" if len(items) > k else ""
                return f"({len(items)}) {head}{tail}"

            print(f"    obs columns : {_short(adata.obs.columns)}")
            print(f"    var columns : {_short(adata.var.columns)}")
            obsm_keys = list(adata.obsm.keys()) if hasattr(adata, "obsm") else []
            print(f"    obsm keys   : {obsm_keys}")
            uns_keys = list(adata.uns.keys()) if hasattr(adata, "uns") else []
            print(f"    uns keys    : {uns_keys}")
            adt_obsm = [k for k in obsm_keys
                        if any(s in k.lower() for s in ("adt", "antibody", "protein"))]
            for k in adt_obsm:
                m = adata.obsm[k]
                print(f"    obsm['{k}'] shape={getattr(m, 'shape', None)}"
                      + (f" columns e.g. {list(m.columns)[:8]}"
                         if hasattr(m, "columns") else ""))
            if "feature_types" in list(adata.var.columns):
                print(f"    var feature_types: "
                      f"{adata.var['feature_types'].value_counts().to_dict()}")

        # --- ENSG list aligned to var_names ---
        ensg_col = find_ensg_column(list(adata.var.columns))
        if ensg_col:
            print(f"    ENSG source : var['{ensg_col}']")
            var_ensg_list = [extract_ensg(v) for v in adata.var[ensg_col].tolist()]
        else:
            print("    ENSG source : var_name strings only (no ENSG var column)")
            var_ensg_list = [extract_ensg(n) for n in var_names]

        # --- Feature-ID list aligned to var_names (grounds custom HTO names) ---
        id_col = find_id_column(list(adata.var.columns))
        if id_col:
            print(f"    ID source     : var['{id_col}']")
            var_ids_list = [
                None if (v is None or (isinstance(v, float) and np.isnan(v))
                         or str(v).strip().lower() == "nan") else str(v)
                for v in adata.var[id_col].tolist()]
        else:
            print("    ID source     : none (no ID var column)")
            var_ids_list = None

        # --- In-matrix feature types (authoritative for THIS X) ---
        var_map = discover_var_feature_types(
            var_names, adata.var["feature_types"]
            if "feature_types" in list(adata.var.columns) else None)
        merged_map, merged_prov = ref_map, ref_provenance
        if var_map:
            merged_map = {**(ref_map or {}), **var_map}  # var column wins
            merged_prov = ("; ".join(p for p in
                                     [ref_provenance,
                                      f"var['feature_types'] ({len(var_map)} typed rows)"]
                                     if p))
        n_nongex = sum(1 for t in (merged_map or {}).values()
                       if t in NON_GEX_FEATURE_TYPES)

        # --- Demux bookkeeping traces (obs/uns) ---
        bookkeeping = summarize_mux_bookkeeping(adata)
        if bookkeeping["class_column"]:
            print(f"    Demux class column: obs['{bookkeeping['class_column']}'] "
                  f"({len(bookkeeping['class_values'])} labels)")
        if bookkeeping["uns_hto_lists"]:
            print(f"    uns tag lists: "
                  f"{ {k: v for k, v in bookkeeping['uns_hto_lists'].items()} }")
        if bookkeeping["uns_filtered_lists"]:
            print(f"    uns filtered lists: "
                  f"{ {k: v for k, v in bookkeeping['uns_filtered_lists'].items()} }")
        if bookkeeping["tag_obs_var_overlap"]:
            print(f"    obs/var tag overlap: {bookkeeping['tag_obs_var_overlap']}")

        # --- Evidence-tier classification ---
        print("\n  [Tiers] Evidence-based classification (standard table)")
        class_df = classify_features(var_names, var_ensg_list,
                                     gencode_symbols, gencode_ensg,
                                     hgnc_cache, ncbi_cache, gmaps,
                                     ref_map=merged_map)
        if merged_map:
            print(f"    Feature reference: {merged_prov} "
                  f"({len(merged_map)} rows, {n_nongex} non-GEX)")
        else:
            print(f"    Feature reference: none discovered ({merged_prov})")
        table = build_standard_table(class_df, ref_provenance=merged_prov)
        if verbose:
            with pd.option_context("display.max_rows", None,
                                   "display.max_columns", None,
                                   "display.width", 250,
                                   "display.max_colwidth", 90):
                print(table.to_string(index=False))
        else:
            core = table[table["category"].str.startswith("match ")
                         | (table["category"] == "unresolved · kept")
                         | table["category"].str.startswith("reference: ")]
            with pd.option_context("display.max_rows", None,
                                   "display.max_columns", None,
                                   "display.width", 150,
                                   "display.max_colwidth", 45):
                print(core.to_string(index=False))
            n_hidden = len(table) - len(core)
            if n_hidden:
                print(f"    ... +{n_hidden} detail rows hidden "
                      f"(biotype/tag breakdown — use --verbose or --out CSV)")

        # --- Distribution stats (optional, read-only) ---
        # Confirmed transcripts (authoritative gene identity) are never
        # HTO/ADT suspects — except corroborated distribution evidence
        # (T-DIST-01), which is reported for awareness, never concealed.
        confirmed_names = set(
            class_df[class_df["confirmed"]]["var_name"].tolist())
        print(f"    Confirmed transcripts: {len(confirmed_names)}/{len(class_df)} "
              f"(exempt from suspicion verdicts; distribution evidence on "
              f"them is still reported)")
        if run_distro:
            print("\n  [Distro] Per-feature counts (read-only)")
            stats_df = compute_distribution_stats(adata.X, var_names,
                                                  chunk_size=chunk_size)
            print("\n  [Uns] Unsupervised hashtag discovery (counts only; "
                  "all features, candidates annotated by identity)")
            unsupervised = discover_hashtag_candidates(stats_df, adata.X,
                                                       var_names)
            if unsupervised["decision"] == "panel":
                print(f"    points at: {unsupervised['panel']}")
            print(f"    {unsupervised['reason']}")
        else:
            print("\n  [Distro] skipped (--no-distro)")
            stats_df = None
            unsupervised = None

        # --- REF tests (run always; need no counts) ---
        print("\n  [Tests:REF] names/symbols/ncRNA vs references")
        ref_tests = run_ref_tests(var_names, class_df, var_ensg_list,
                                  ensg_col, gencode_ensg, gmaps)

        # --- DIST tests (need counts; minimal fallback covers ADT-like names) ---
        adt_like_names = [n for n in var_names if "adt-like" in tag_name(n)[1]]
        hto_like_names = [n for n in var_names if "hto-like" in tag_name(n)[1]]
        # Feature-barcode names from the merged reference (file + var column)
        # must have counts available even in minimal mode (T-REF-09 needs them).
        ref_bc_in_var = ([n for n in (merged_map or {})
                          if merged_map[n] in NON_GEX_FEATURE_TYPES and n in var_names]
                         if merged_map else [])
        # HTO-ID names likewise (catches HTOs regardless of their ref type).
        hto_id_in_var = [n for n, i in zip(var_names, var_ids_list or [])
                         if i and HTO_ID_RE.search(str(i))]
        if stats_df is None:
            mini = []
            for n in sorted(set(adt_like_names) | set(hto_like_names)
                            | set(ref_bc_in_var) | set(hto_id_in_var)):
                j = var_names.index(n)
                col = _dense_slice(adata.X, j, j + 1).flatten()
                mini.append({"var_name": n, "mean": float(col.mean()),
                             "var": 0.0, "cv": 0.0, "overdisp": 0.0,
                             "frac_expressed": float((col > 0).mean()),
                             "max": float(col.max()), "bc": float("nan")})
            stats_df = pd.DataFrame(
                mini, columns=["var_name", "mean", "var", "cv", "overdisp",
                               "frac_expressed", "max", "bc"])
            print("\n  [Tests:DIST] minimal feature-barcode counts (full distro skipped)")
            dist_tests, verdict = run_dist_tests(
                class_df, stats_df, adt_like_names, adt_panel,
                full_distro=False, confirmed_names=confirmed_names)
        else:
            print("\n  [Tests:DIST] counts vs empirical ADT panel")
            panel_names = (set(unsupervised["panel"])
                           if unsupervised and unsupervised.get("decision") == "panel"
                           else set())
            dist_tests, verdict = run_dist_tests(class_df, stats_df,
                                                 adt_like_names, adt_panel,
                                                 full_distro=True,
                                                 confirmed_names=confirmed_names,
                                                 panel_names=panel_names)

        # --- MUX tests (HTO/multiplexing tags + reference reconciliation) ---
        print("\n  [Tests:MUX] multiplexing tags vs feature reference")
        obsm_col_sets = {}
        for k in obsm_keys:
            m = adata.obsm[k]
            if hasattr(m, "columns"):
                obsm_col_sets[k] = list(m.columns)
        mux_tests, mux_placement = run_mux_tests(
            var_names, class_df, stats_df, obsm_col_sets, merged_map,
            merged_prov, full_distro=run_distro, bookkeeping=bookkeeping,
            var_ids_list=var_ids_list, unsupervised=unsupervised)

        all_tests = ref_tests + dist_tests + mux_tests
        test_df = pd.DataFrame(
            [{"test": t, "check": test_title(t), "domain": d, "result": r,
              "summary": s, "n_affected": n,
              "examples": "; ".join(map(str, ex[:5]))}
             for t, d, r, s, n, ex in all_tests],
            columns=["test", "check", "domain", "result", "summary",
                     "n_affected", "examples"])
        if verbose:
            with pd.option_context("display.max_rows", None,
                                   "display.max_columns", None,
                                   "display.width", 250,
                                   "display.max_colwidth", 120):
                print(test_df.to_string(index=False))
        else:
            n_pass = sum(1 for t in all_tests if t[2] == "PASS")
            print(f"    ... {n_pass}/{len(all_tests)} checks passed "
                  f"(only WARN/FAIL shown below)")
            for t, d, r, s, n, ex in all_tests:
                if r == "PASS":
                    continue
                short = s if len(s) <= 220 else s[:217] + "..."
                print(f"    [{r}] {test_title(t)} ({t}) — {short}")
        # --- Explicit ADT placement declaration: where antibody signals live ---
        adt_hits = (stats_df[stats_df["var_name"].isin(adt_like_names)]
                    if stats_df is not None and adt_like_names
                    else pd.DataFrame(columns=["mean"]))
        n_adt_nonzero = int((adt_hits["mean"] > 0).sum()) if len(adt_hits) else 0
        obsm_loc = ", ".join(f"obsm['{k}']" for k in adt_obsm) or "—"
        if adt_obsm and not adt_like_names:
            placement = (f"ADTs are present and are placed in {obsm_loc}; "
                         "absent from X/var (0 ADT-like var_names).")
        elif adt_obsm and adt_like_names:
            placement = (f"ADTs are present and are placed in {obsm_loc}; "
                         f"additionally {len(adt_like_names)} ADT-like var_names in "
                         f"X/var ({n_adt_nonzero} nonzero, "
                         f"{len(adt_like_names) - n_adt_nonzero} zero-count).")
        elif not adt_obsm and adt_like_names:
            placement = (f"ADTs are present and are placed in X/var "
                         f"({len(adt_like_names)} ADT-like var_names, "
                         f"{n_adt_nonzero} nonzero, "
                         f"{len(adt_like_names) - n_adt_nonzero} zero-count placeholders); "
                         "no antibody matrix in obsm.")
        else:
            placement = ("No antibody signals detected: no ADT matrix in obsm "
                         "and no ADT-like var_names in X/var.")
        print = _print  # verdicts always print, every mode
        if brief:
            # Ultra-brief: just the essential verdict and key findings
            fails = test_df[test_df["result"] == "FAIL"]["test"].tolist()
            warns = test_df[test_df["result"] == "WARN"]["test"].tolist()
            fail_names = ",".join(test_title(t) for t in fails)
            warn_names = ",".join(test_title(t) for t in warns)
            if fails:
                overall = (f"FAIL ({fail_names}) — issue found, "
                           "curator review required, nothing was removed")
            elif warns:
                overall = (f"PASS WITH NOTES ({warn_names}) — "
                           "no corruption found, review lists above are retained data")
            else:
                overall = "CLEAN — all tests passed"
            print(f"  OVERALL: {overall}")
            print(f"  ADT: {placement}")
            print(f"  MUX: {mux_placement}")
            # Only show WARN/FAIL one-liners
            for t, d, r, s, n, ex in all_tests:
                if r in ("WARN", "FAIL"):
                    short = s if len(s) <= 160 else s[:157] + "..."
                    print(f"  [{r}] {test_title(t)} — {short}")
        else:
            print(f"\n    ADT PLACEMENT: {placement}")
            print(f"\n    ADT VERDICT: {verdict}")
            print(f"\n    MUX PLACEMENT: {mux_placement}")
            # --- Plain-language overall verdict (the one line to read) ---
            fails = test_df[test_df["result"] == "FAIL"]["test"].tolist()
            warns = test_df[test_df["result"] == "WARN"]["test"].tolist()
            fail_names = ",".join(test_title(t) for t in fails)
            warn_names = ",".join(test_title(t) for t in warns)
            if fails:
                overall = (f"FAIL ({fail_names}) — issue found, "
                           "curator review required, nothing was removed")
            elif warns:
                overall = (f"PASS WITH NOTES ({warn_names}) — "
                           "no corruption found, review lists above are retained data")
            else:
                overall = "CLEAN — all tests passed"
            print(f"\n    OVERALL: {overall}")
    finally:
        adata.file.close()

    return {
        "dataset": label,
        "h5ad_path": h5ad_path,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "classification": class_df,
        "standard_table": table,
        "test_table": test_df,
        "adt_placement": placement,
        "mux_placement": mux_placement,
        "adt_verdict": verdict,
        "overall": overall,
    }


# ============================================================================
# HGNC online refresh (opt-in only; never runs by default)
# ============================================================================

def refresh_hgnc_cache(symbols, cache_dir=None):
    """Query HGNC REST for symbols missing from the cache. Opt-in only."""
    import requests
    if cache_dir is None:
        cache_dir = CACHE_DIR
    cache = load_hgnc_cache(cache_dir)
    missing = [s for s in symbols if s not in cache]
    print(f"HGNC refresh: {len(missing)} symbols missing from cache")
    for i, sym in enumerate(missing):
        try:
            resp = requests.get(f"{HGNC_API}/search/{sym}",
                                headers={"Accept": "application/json"},
                                timeout=15)
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
            with open(os.path.join(cache_dir, "hgnc_mapping.json"), "w") as f:
                json.dump(cache, f, indent=2)
            print(f"  ... {i + 1}/{len(missing)}")
    with open(os.path.join(cache_dir, "hgnc_mapping.json"), "w") as f:
        json.dump(cache, f, indent=2)
    print("HGNC cache saved.")
