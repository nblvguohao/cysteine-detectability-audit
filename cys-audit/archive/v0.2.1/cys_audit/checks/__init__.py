"""The five artefact tests. Each `run(ds, cfg)` returns a dict with at least
status, reason, statistic, null, margin, estimate, ci, n_positive, n_background, and details."""

from . import abundance, cleavage, multicys, overlap, searchspace  # noqa: F401

ORDER = ("cleavage_geometry", "multi_cysteine", "protein_abundance", "positive_detected_overlap", "search_space")
