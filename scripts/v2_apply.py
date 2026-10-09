"""Apply v2 models to sequences outside the training matrix.

Two consumers use this module: the retrospective external audit (human
PXD044043, rice PXD072089) and a candidate table.  Neither cohort has
AlphaFold-derived structural features available here, so the five structural
columns of the frozen block are passed as missing, exactly as in the v1
external audit.  ProtBERT embeddings were produced on a separate CUDA host and
are not available for these sequences, so the deployed model is the single
LightGBM binary member rather than the three-member blend.
"""
from __future__ import annotations

import importlib.util
import json
import os
import statistics
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from common import INPUTS, RESULTS, ROOT
import v2_features as features
import v2_members as members
from v2_stack import load_v2_matrix

DEFAULT_SOURCE_PROJECT = Path(
    "/path/to/source_project"
)
SOURCE_PROJECT = Path(
    os.environ.get("PERSULFIDATION_SOURCE_PROJECT", DEFAULT_SOURCE_PROJECT)
)
SOURCE_PIPELINE = SOURCE_PROJECT / "outputs/20260911_模型与Q1c推进/model_pipeline.py"
DETECTABILITY_TOKENS = ("T:", "nearest_K", "nearest_R", "pep_", "kr_within")


def load_source_encoder():
    spec = importlib.util.spec_from_file_location("source_model_pipeline", SOURCE_PIPELINE)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"Cannot load {SOURCE_PIPELINE}")
    sys.path.insert(0, str(SOURCE_PIPELINE.parent))
    spec.loader.exec_module(module)
    if len(module.feature_names("all_basic_removed")) != 23:
        raise RuntimeError("Unexpected frozen feature width from the source encoder")
    return module


def build_matrix(records, encoder, names_reference, structure_lookup=None):
    """Featurise new sequences in the v2 layout.

    ``records`` is a sequence of dicts with ``accession`` and ``sequence``.
    Every cysteine of every sequence becomes one row.  ``structure_lookup``
    optionally supplies the per-site structural dictionary that the source
    encoder expects; when it is absent the five structural columns are missing,
    exactly as in the v1 external audit.
    """
    proteins = []
    positions = []
    base_rows = []
    site_names = None
    protein_names = None
    for record in records:
        sequence = features._normalise(record["sequence"].strip().upper())
        arrays = features._protein_arrays(sequence)
        protein_block_names, protein_values = features._protein_features(sequence)
        for index, residue in enumerate(sequence):
            if residue != "C":
                continue
            position = index + 1
            structure = None
            if structure_lookup is not None:
                structure = structure_lookup.get(record["accession"], {}).get(position)
            frozen = encoder.vector(sequence, position, structure, "all_basic_removed")
            new_names, new_values = features.site_base_block(sequence, position, arrays)
            if site_names is None:
                site_names = new_names
                protein_names = protein_block_names
            elif new_names != site_names or protein_block_names != protein_names:
                raise RuntimeError("Feature name order changed between sites")
            base_rows.append(list(frozen) + list(new_values) + list(protein_values))
            proteins.append(record["accession"])
            positions.append(position)
    proteins = np.asarray(proteins)
    positions = np.asarray(positions, dtype=int)
    base = np.asarray(base_rows, dtype=np.float32)
    base_names = [f"A:frozen_{i:02d}" for i in range(23)] + list(site_names) + list(protein_names)
    context, context_names = features.context_block(base, base_names, proteins)
    matrix = np.hstack([base, context])
    matrix_names = base_names + context_names
    if matrix_names != list(names_reference):
        raise RuntimeError("External feature layout does not match the training layout")
    return matrix, proteins, positions


def select_columns(names, feature_set):
    if feature_set == "full":
        return np.arange(len(names))
    if feature_set == "chem":
        return np.asarray([
            i for i, name in enumerate(names)
            if not any(token in name for token in DETECTABILITY_TOKENS)
        ])
    raise SystemExit(f"Unknown feature set {feature_set}")


def deployed_configuration(feature_set):
    """Median boosting rounds and majority configuration across outer folds.

    Taken from the stack run's per-fold selections, which were made on inner
    validation folds only.  No external label is consulted.
    """
    path = RESULTS / f"v2_stack_{feature_set}_fold_records.json"
    records = json.loads(path.read_text())
    rounds = []
    configs = []
    for record in records:
        selection = record["members"]["lgb_bin"]["selection"]
        rounds.append(int(selection["rounds"]))
        configs.append(json.dumps(selection["config"], sort_keys=True))
    majority = Counter(configs).most_common(1)[0][0]
    return {
        "config": json.loads(majority),
        "rounds": int(statistics.median(rounds)),
        "per_fold_rounds": rounds,
        "per_fold_configs": [json.loads(config) for config in configs],
        "source": str(path.relative_to(ROOT)),
    }


def fit_deployed_model(feature_set, seeds=None):
    """Fit the deployed LightGBM binary member on every site of the training matrix.

    With ``seeds`` a small ensemble is fitted, differing only in the LightGBM
    seeds, so that the spread of a site score across seeds can be reported as a
    stability figure.  The boosting-round count and configuration always come
    from the inner-fold selections of the stack run.
    """
    X, names, y, folds, proteins, components, positions, meta = load_v2_matrix()
    columns = select_columns(names, feature_set)
    Xc = np.ascontiguousarray(X[:, columns])
    deployment = deployed_configuration(feature_set)
    boosters = []
    for seed in (seeds or [None]):
        config = dict(deployment["config"])
        params = members._params("bin", config)
        if seed is not None:
            params.update({
                "seed": int(seed),
                "feature_fraction_seed": int(seed) + 1,
                "bagging_seed": int(seed) + 2,
                "data_random_seed": int(seed) + 3,
            })
        data, _ = members._dataset(Xc, y, np.arange(len(y)), proteins, "bin")
        boosters.append(
            members.lgb.train(params, data, num_boost_round=deployment["rounds"])
        )
    if seeds is None:
        return boosters[0], names, columns, deployment
    return boosters, names, columns, deployment
