"""scverify — read-only verification that ADT/HTO counts are not mixed into
GEX matrices of single-cell .h5ad files.

Public API:
  verify_dataset(...)            full verification of one dataset
  classify_features(...)         evidence-tier classification of var_names
  compute_distribution_stats(...) per-feature count statistics
  discover_feature_reference(...) authoritative feature types from raw/ files
  discover_var_feature_types(...) feature types from a var column
  summarize_mux_bookkeeping(...) demultiplexing traces in obs/uns
  discover_hashtag_candidates(...) unsupervised hashtag search
  load_all_caches(...)           every authority cache in one dict
  refresh_hgnc_cache(...)        opt-in HGNC REST refresh of the local cache
"""

from scverify.verification import (  # noqa: F401
    CACHE_DIR,
    HGNC_API,
    NCRNA_BIOTYPES,
    PATTERNS,
    classify_features,
    compute_distribution_stats,
    discover_feature_reference,
    discover_hashtag_candidates,
    discover_var_feature_types,
    load_adt_panel,
    load_all_caches,
    load_gencode_ensg,
    load_gencode_maps,
    load_gencode_symbols,
    load_hgnc_cache,
    load_ncbi_cache,
    refresh_hgnc_cache,
    summarize_mux_bookkeeping,
    tag_name,
    test_title,
    verify_dataset,
)

__version__ = "0.1.0"
