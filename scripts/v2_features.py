"""Label-free v2 feature construction for Cys annotation ranking.

Three blocks are produced in the frozen site order of
``inputs/primary_site_folds.csv``:

* ``A`` the unchanged 23 frozen engineered features.
* ``B`` new multi-scale sequence-context features derived only from the tomato
  proteome sequence and the Cys position.
* ``C`` within-protein contextual transforms (rank percentile and protein-mean
  centring) of a predeclared subset of ``A`` and ``B``.

No label, fold index or evaluation metric is read anywhere in this module, so
the transform is identical inside and outside any training partition.  Every
protein belongs to exactly one homology component and therefore to exactly one
outer fold, so within-protein statistics cannot move information across folds.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

from common import INPUTS, ROOT, read_csv, load_frozen, sha256, write_json

FEATURES = ROOT / "features"

AMINO = "ACDEFGHIKLMNPQRSTVWY"

# Kyte & Doolittle 1982 hydropathy.
HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
    "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
    "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
# Formal side-chain charge near pH 7.
CHARGE = {a: 0.0 for a in AMINO}
CHARGE.update({"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0})
# Zamyatnin 1972 residue volumes (A^3).
VOLUME = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1,
    "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1,
    "P": 112.7, "Q": 143.8, "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0,
    "W": 227.8, "Y": 193.6,
}
# Bhaskaran & Ponnuswamy 1988 average flexibility.
FLEXIBILITY = {
    "A": 0.357, "C": 0.346, "D": 0.511, "E": 0.497, "F": 0.314, "G": 0.544,
    "H": 0.323, "I": 0.462, "K": 0.466, "L": 0.365, "M": 0.295, "N": 0.463,
    "P": 0.509, "Q": 0.493, "R": 0.529, "S": 0.507, "T": 0.444, "V": 0.386,
    "W": 0.305, "Y": 0.420,
}
# Grantham 1974 polarity.
POLARITY = {
    "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3, "F": 5.2, "G": 9.0, "H": 10.4,
    "I": 5.2, "K": 11.3, "L": 4.9, "M": 5.7, "N": 11.6, "P": 8.0, "Q": 10.5,
    "R": 10.5, "S": 9.2, "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2,
}
# Chou & Fasman 1978 conformational parameters.
HELIX = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13, "G": 0.57,
    "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21, "M": 1.45, "N": 0.67,
    "P": 0.57, "Q": 1.11, "R": 0.98, "S": 0.77, "T": 0.83, "V": 1.06,
    "W": 1.08, "Y": 0.69,
}
SHEET = {
    "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38, "G": 0.75,
    "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30, "M": 1.05, "N": 0.89,
    "P": 0.55, "Q": 1.10, "R": 0.93, "S": 0.75, "T": 1.19, "V": 1.70,
    "W": 1.37, "Y": 1.47,
}

SCALES = {
    "hydropathy": HYDROPATHY,
    "charge": CHARGE,
    "volume": VOLUME,
    "flexibility": FLEXIBILITY,
    "polarity": POLARITY,
    "helix": HELIX,
    "sheet": SHEET,
}

CLASSES = {
    "hydrophobic": set("AVILMFWC"),
    "aromatic": set("FWYH"),
    "acidic": set("DE"),
    "basic": set("KRH"),
    "polar": set("STNQCYH"),
    "small": set("AGSCTDNP"),
    "tiny": set("AGS"),
    "gly": set("G"),
    "pro": set("P"),
    "cys": set("C"),
    "serthr": set("ST"),
    "his": set("H"),
    "met": set("M"),
    "aliphatic": set("AVILM"),
}

WINDOWS = (1, 2, 3, 5, 7, 10, 15, 20, 30, 50)
OFFSETS = tuple(o for o in range(-8, 9) if o != 0)
OFFSET_SCALES = ("hydropathy", "charge", "volume", "flexibility", "polarity")
INDICATOR_OFFSETS = tuple(o for o in range(-4, 5) if o != 0)
INDICATOR_CLASSES = ("basic", "acidic", "aromatic", "gly", "pro", "serthr", "cys")
CYS_PARTNER_OFFSETS = tuple(o for o in range(-6, 7) if o != 0)

# Predeclared subset receiving within-protein contextual transforms.  The rule
# is: every frozen feature that varies between Cys of one protein, plus the
# multi-scale window and Cys-topology features that the within-protein ranking
# metric depends on.  Fixed before any v2 metric was computed.
CONTEXT_PREFIXES = (
    "A:",
    "T:pep_len", "T:pep_mass", "T:pep_gravy", "T:pep_detectable", "T:pep_mc",
    "T:pep_relative_offset", "T:pep_kr_within_20",
    "B:win5_", "B:win10_", "B:win20_",
    "B:cys_dist_", "B:cys_in_", "B:nearest_",
    "B:pos_",
)
CONTEXT_EXCLUDE_SUBSTRINGS = ("_length", "protein_", "window_len")


def _sequences():
    sequences = {}
    with (INPUTS / "sly_proteome.tsv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequences[row["Entry"]] = row["Sequence"].strip().upper()
    return sequences


def _normalise(sequence):
    return "".join(residue if residue in AMINO else "X" for residue in sequence)


def _protein_arrays(sequence):
    """Per-protein residue arrays reused by every Cys site of that protein."""
    scale_arrays = {
        name: np.asarray([table.get(residue, 0.0) for residue in sequence], dtype=np.float64)
        for name, table in SCALES.items()
    }
    class_arrays = {
        name: np.asarray([residue in members for residue in sequence], dtype=np.float64)
        for name, members in CLASSES.items()
    }
    letters = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    boundaries = _cleavage_boundaries(sequence)
    residue_positions = {
        residue: np.flatnonzero(letters == residue.encode("ascii"))
        for residue in ("C", "H", "M", "W", "K", "R", "D", "E")
    }
    return {
        "scales": scale_arrays,
        "classes": class_arrays,
        "residue_positions": residue_positions,
        "boundaries": boundaries,
    }


def _site_features(sequence, position, arrays):
    """Return (names, values) for one Cys site; ``position`` is 1-based."""
    names = []
    values = []
    n = len(sequence)
    index = position - 1
    scale_arrays = arrays["scales"]
    class_arrays = arrays["classes"]
    residue_positions = arrays["residue_positions"]

    for width in WINDOWS:
        low = max(0, index - width)
        high = min(n, index + width + 1)
        flank = list(range(low, index)) + list(range(index + 1, high))
        flank_index = np.asarray(flank, dtype=int)
        size = len(flank_index)
        prefix = f"B:win{width}_"
        for scale_name in SCALES:
            names.append(prefix + f"mean_{scale_name}")
            values.append(
                float(scale_arrays[scale_name][flank_index].mean()) if size else 0.0
            )
        for class_name in CLASSES:
            names.append(prefix + f"frac_{class_name}")
            values.append(
                float(class_arrays[class_name][flank_index].mean()) if size else 0.0
            )
        for class_name in ("cys", "his", "acidic", "basic"):
            names.append(prefix + f"count_{class_name}")
            values.append(
                float(class_arrays[class_name][flank_index].sum()) if size else 0.0
            )
        names.append(prefix + "net_charge")
        values.append(float(scale_arrays["charge"][flank_index].sum()) if size else 0.0)
        names.append(prefix + "window_len")
        values.append(float(size))

    for offset in OFFSETS:
        target = index + offset
        inside = 0 <= target < n
        for scale_name in OFFSET_SCALES:
            names.append(f"B:pos_{offset:+d}_{scale_name}")
            values.append(float(scale_arrays[scale_name][target]) if inside else 0.0)
        names.append(f"B:pos_{offset:+d}_present")
        values.append(1.0 if inside else 0.0)

    for offset in INDICATOR_OFFSETS:
        target = index + offset
        inside = 0 <= target < n
        for class_name in INDICATOR_CLASSES:
            names.append(f"B:pos_{offset:+d}_is_{class_name}")
            values.append(float(class_arrays[class_name][target]) if inside else 0.0)

    cys_positions = residue_positions["C"]
    other = cys_positions[cys_positions != index]
    left = other[other < index]
    right = other[other > index]
    names.append("B:cys_dist_left")
    values.append(float(index - left.max()) if left.size else float(n))
    names.append("B:cys_dist_right")
    values.append(float(right.min() - index) if right.size else float(n))
    names.append("B:cys_dist_min")
    values.append(
        float(np.abs(other - index).min()) if other.size else float(n)
    )
    for radius in (4, 8, 16, 32, 64):
        names.append(f"B:cys_in_{radius}")
        values.append(float(np.sum(np.abs(other - index) <= radius)))
    for offset in CYS_PARTNER_OFFSETS:
        target = index + offset
        names.append(f"B:cys_partner_{offset:+d}")
        values.append(1.0 if 0 <= target < n and sequence[target] == "C" else 0.0)

    for residue in ("H", "M", "W", "K", "R", "D", "E"):
        hits = residue_positions[residue]
        hits = hits[hits != index]
        names.append(f"B:nearest_{residue}")
        values.append(float(np.abs(hits - index).min()) if hits.size else float(n))

    names.append("B:pos_index_1based")
    values.append(float(position))
    names.append("B:pos_from_n_term")
    values.append(float(index))
    names.append("B:pos_from_c_term")
    values.append(float(n - 1 - index))
    names.append("B:pos_min_terminus")
    values.append(float(min(index, n - 1 - index)))
    names.append("B:pos_relative")
    values.append(float(index) / float(max(1, n - 1)))
    names.append("B:cys_rank_in_protein")
    values.append(float(np.searchsorted(cys_positions, index) + 1))
    names.append("B:cys_rank_fraction")
    values.append(
        float(np.searchsorted(cys_positions, index) + 1) / float(max(1, cys_positions.size))
    )
    return names, values


# Monoisotopic residue masses (Da) for in-silico tryptic peptide descriptors.
RESIDUE_MASS = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
    "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
    "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
    "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
    "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
    "X": 110.0,
}
WATER_MASS = 18.010565
DETECTABLE_LENGTH = (7, 30)
DETECTABLE_MASS = (700.0, 3500.0)


def _cleavage_boundaries(sequence):
    """Trypsin boundaries: cut after K/R unless followed by P."""
    n = len(sequence)
    boundaries = [0]
    for i, residue in enumerate(sequence):
        if residue in ("K", "R") and (i + 1 >= n or sequence[i + 1] != "P"):
            boundaries.append(i + 1)
    if boundaries[-1] != n:
        boundaries.append(n)
    return boundaries


def _peptide_descriptors(sequence, start, stop):
    peptide = sequence[start:stop]
    length = len(peptide)
    if length == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    mass = WATER_MASS + sum(RESIDUE_MASS.get(r, RESIDUE_MASS["X"]) for r in peptide)
    gravy = float(np.mean([HYDROPATHY.get(r, 0.0) for r in peptide]))
    charge = float(sum(CHARGE.get(r, 0.0) for r in peptide))
    cys = float(peptide.count("C"))
    return float(length), float(mass), gravy, charge, cys


def _digest_features(sequence, position, boundaries):
    """In-silico tryptic peptide descriptors for one Cys site.

    These describe how detectable the Cys-containing peptide is in a
    trypsin-based proteomics workflow.  They are deliberately separated into
    their own block so that the dependence of the annotation label on method
    detectability can be measured by ablation instead of being hidden inside
    the sequence block.
    """
    import bisect

    names = []
    values = []
    n = len(sequence)
    index = position - 1
    peptide_index = bisect.bisect_right(boundaries, index) - 1
    peptide_index = min(max(peptide_index, 0), len(boundaries) - 2)
    start = boundaries[peptide_index]
    stop = boundaries[peptide_index + 1]

    length, mass, gravy, charge, cys = _peptide_descriptors(sequence, start, stop)
    names += [
        "T:pep_len", "T:pep_mass", "T:pep_gravy", "T:pep_charge", "T:pep_cys_count",
    ]
    values += [length, mass, gravy, charge, cys]
    names.append("T:pep_offset_in_peptide")
    values.append(float(index - start))
    names.append("T:pep_relative_offset")
    values.append(float(index - start) / float(max(1.0, length - 1)))
    names.append("T:pep_is_protein_n_terminal")
    values.append(1.0 if start == 0 else 0.0)
    names.append("T:pep_is_protein_c_terminal")
    values.append(1.0 if stop == n else 0.0)
    names.append("T:pep_detectable_length")
    values.append(
        1.0 if DETECTABLE_LENGTH[0] <= length <= DETECTABLE_LENGTH[1] else 0.0
    )
    names.append("T:pep_detectable_mass")
    values.append(1.0 if DETECTABLE_MASS[0] <= mass <= DETECTABLE_MASS[1] else 0.0)
    names.append("T:pep_log_len")
    values.append(float(np.log1p(length)))

    detectable_any = values[names.index("T:pep_detectable_length")]
    for missed in (1, 2):
        left = boundaries[max(0, peptide_index - missed)]
        right = boundaries[min(len(boundaries) - 1, peptide_index + 1 + missed)]
        left_length, left_mass, _, _, _ = _peptide_descriptors(sequence, left, stop)
        right_length, right_mass, _, _, _ = _peptide_descriptors(sequence, start, right)
        names.append(f"T:pep_mc{missed}_left_len")
        values.append(left_length)
        names.append(f"T:pep_mc{missed}_right_len")
        values.append(right_length)
        names.append(f"T:pep_mc{missed}_min_len")
        values.append(min(left_length, right_length))
        detectable = 0.0
        for candidate_length, candidate_mass in (
            (left_length, left_mass), (right_length, right_mass)
        ):
            if (
                DETECTABLE_LENGTH[0] <= candidate_length <= DETECTABLE_LENGTH[1]
                and DETECTABLE_MASS[0] <= candidate_mass <= DETECTABLE_MASS[1]
            ):
                detectable = 1.0
        names.append(f"T:pep_mc{missed}_detectable")
        values.append(detectable)
        detectable_any = max(detectable_any, detectable)
    names.append("T:pep_detectable_any_missed_cleavage")
    values.append(detectable_any)

    names.append("T:pep_n_peptides_in_protein")
    values.append(float(len(boundaries) - 1))
    names.append("T:pep_kr_within_20")
    values.append(
        float(sum(1 for r in sequence[max(0, index - 20):min(n, index + 21)] if r in "KR"))
    )
    names.append("T:pep_has_kp_or_rp")
    values.append(
        1.0 if any(
            sequence[i] in "KR" and i + 1 < n and sequence[i + 1] == "P"
            for i in range(start, stop)
        ) else 0.0
    )
    names.append("T:pep_met_count")
    values.append(float(sequence[start:stop].count("M")))
    return names, values


def _protein_features(sequence):
    names = []
    values = []
    n = len(sequence)
    array = np.frombuffer(sequence.encode("ascii"), dtype="S1")
    names.append("B:protein_length")
    values.append(float(n))
    names.append("B:protein_log_length")
    values.append(float(np.log1p(n)))
    for residue in AMINO:
        names.append(f"B:protein_frac_{residue}")
        values.append(float(np.mean(array == residue.encode("ascii"))))
    hydropathy = np.asarray([HYDROPATHY.get(r, 0.0) for r in sequence])
    charge = np.asarray([CHARGE.get(r, 0.0) for r in sequence])
    names.append("B:protein_mean_hydropathy")
    values.append(float(hydropathy.mean()))
    names.append("B:protein_net_charge")
    values.append(float(charge.sum()))
    names.append("B:protein_charge_density")
    values.append(float(charge.sum()) / float(n))
    cys = int(np.sum(array == b"C"))
    names.append("B:protein_cys_count")
    values.append(float(cys))
    names.append("B:protein_cys_density")
    values.append(float(cys) / float(n))
    return names, values


def context_columns(base_names):
    """Predeclared columns that receive within-protein contextual transforms."""
    keep = []
    for j, name in enumerate(base_names):
        if not name.startswith(CONTEXT_PREFIXES):
            continue
        if any(token in name for token in CONTEXT_EXCLUDE_SUBSTRINGS):
            continue
        keep.append(j)
    return np.asarray(keep, dtype=int)


def context_block(base, base_names, group_labels):
    """Within-protein rank percentile and protein-mean centring.

    Uses feature values only.  Every protein sits inside a single homology
    component and therefore inside a single outer fold, so these statistics
    cannot carry information between folds.
    """
    keep = context_columns(base_names)
    half = len(keep)
    context = np.zeros((base.shape[0], 2 * half), dtype=np.float32)
    context_names = [f"C:rank_{base_names[j]}" for j in keep] + [
        f"C:centred_{base_names[j]}" for j in keep
    ]
    order = np.argsort(group_labels, kind="stable")
    starts = np.flatnonzero(
        np.r_[True, group_labels[order][1:] != group_labels[order][:-1]]
    )
    for index_group in np.split(order, starts[1:]):
        block = base[np.ix_(index_group, keep)].astype(np.float64)
        if block.shape[0] == 1:
            context[index_group, :half] = 0.5
            context[index_group, half:] = 0.0
            continue
        finite = np.isfinite(block)
        ranks = np.full(block.shape, 0.5, dtype=np.float64)
        for column in range(block.shape[1]):
            values = block[:, column]
            observed = np.isfinite(values)
            if observed.sum() < 2:
                continue
            sub = values[observed]
            unique, inverse, counts = np.unique(
                sub, return_inverse=True, return_counts=True
            )
            cumulative = np.cumsum(counts) - counts
            average = cumulative + (counts - 1) / 2.0
            ranks[observed, column] = average[inverse] / float(sub.shape[0] - 1)
        means = np.nanmean(np.where(finite, block, np.nan), axis=0)
        means = np.where(np.isfinite(means), means, 0.0)
        context[np.ix_(index_group, np.arange(half))] = ranks.astype(np.float32)
        context[np.ix_(index_group, np.arange(half, 2 * half))] = np.where(
            finite, block - means, np.nan
        ).astype(np.float32)
    return context, context_names


def site_base_block(sequence, position, arrays):
    """New sequence-context and digest features for one Cys site."""
    names, values = _site_features(sequence, position, arrays)
    digest_names, digest_values = _digest_features(
        sequence, position, arrays["boundaries"]
    )
    return names + digest_names, values + digest_values


def build():
    started = time.time()
    X_frozen, y, folds, proteins, components, positions, meta = load_frozen()
    sequences = {key: _normalise(value) for key, value in _sequences().items()}
    missing = sorted({p for p in proteins if p not in sequences})
    if missing:
        raise SystemExit(f"Missing sequences for {len(missing)} proteins: {missing[:5]}")

    site_names = None
    protein_names = None
    protein_cache = {}
    collected = []
    for i, (accession, position) in enumerate(zip(proteins, positions)):
        sequence = sequences[accession]
        if sequence[position - 1] != "C":
            raise SystemExit(f"{accession} position {position} is not cysteine")
        if accession not in protein_cache:
            protein_cache[accession] = (
                _protein_features(sequence),
                _protein_arrays(sequence),
            )
        (protein_name_block, protein_values), arrays = protein_cache[accession]
        names, values = site_base_block(sequence, int(position), arrays)
        if site_names is None:
            site_names = names
            protein_names = protein_name_block
        elif names != site_names or protein_name_block != protein_names:
            raise SystemExit("Feature name order changed between sites")
        collected.append(values + protein_values)
        if (i + 1) % 10000 == 0:
            print(f"featurised {i + 1} sites", flush=True)

    B = np.asarray(collected, dtype=np.float32)
    B_names = list(site_names) + list(protein_names)
    A_names = [f"A:frozen_{i:02d}" for i in range(X_frozen.shape[1])]
    base = np.hstack([X_frozen.astype(np.float32), B])
    base_names = A_names + B_names

    context, context_names = context_block(base, base_names, proteins)

    matrix = np.hstack([base, context])
    names = base_names + context_names
    if matrix.shape[1] != len(names):
        raise SystemExit("Name/column mismatch")

    FEATURES.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        FEATURES / "v2_features.npz",
        X=matrix,
        names=np.asarray(names),
        proteins=proteins,
        positions=positions,
        folds=folds,
        components=components,
        y=y,
    )
    audit = {
        "completed": True,
        "elapsed_seconds": time.time() - started,
        "n_sites": int(matrix.shape[0]),
        "n_features_total": int(matrix.shape[1]),
        "n_features_frozen_A": len(A_names),
        "n_features_new_sequence_B": len(B_names),
        "n_features_digest_T": sum(1 for name in B_names if name.startswith("T:")),
        "n_features_within_protein_C": len(context_names),
        "labels_used_in_construction": False,
        "folds_used_in_construction": False,
        "within_protein_statistics": "rank percentile and protein-mean centring over the predeclared subset; every protein lies inside one homology component and one outer fold",
        "context_prefixes": list(CONTEXT_PREFIXES),
        "context_exclude_substrings": list(CONTEXT_EXCLUDE_SUBSTRINGS),
        "windows": list(WINDOWS),
        "offsets": list(OFFSETS),
        "scales": sorted(SCALES),
        "uniprot_keywords_used": False,
        "uniprot_keyword_exclusion_reason": "annotation-derived keywords could encode the label",
        "versions": {"python": sys.version, "numpy": np.__version__},
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                INPUTS / "sly_proteome.tsv",
            ]
        },
        "output_sha256": sha256(FEATURES / "v2_features.npz"),
    }
    write_json(FEATURES / "v2_feature_audit.json", audit)
    print(json.dumps({k: v for k, v in audit.items() if k != "input_hashes"}, indent=2)[:1200])


if __name__ == "__main__":
    build()
