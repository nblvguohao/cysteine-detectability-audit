"""Artefact 5 - search space (chemical identity decided by what the search allowed).

Question: could spectra of an isobaric alternative modification have been assigned to the claimed one?
Two deltas closer than the precursor tolerance cannot be told apart by precursor mass, and if the
alternative was not in the search space it could not be chosen at all.

Built-in table (monoisotopic, from atomic masses in constants.py; extend with care):
  persulfidation (+S)          vs dioxidation / sulfinic acid (+2O)       delta 0.017758 Da
  dioxidation (+2O)            vs persulfidation (+S)                     delta 0.017758 Da
  trioxidation (+3O)           vs oxidised persulfide (+S+O)              delta 0.017758 Da
Resolving mass M* = delta / (tolerance_ppm * 1e-6): below M* the alternative falls outside the
claimed modification's precursor window; at or above M* it falls inside.

Rules (no statistics; configuration plus optional observed peptide masses):
  UNDECIDABLE  identity readout is 'chemistry_inferred' (identity assigned by the enrichment
               chemistry, not by a measured mass - no search space can confirm it); or the
               precursor tolerance is not supplied; or the modification is not in the table
               (reported as 'not covered', never as PASS)
  PASS         alternative in the search space AND M* above the largest peptide mass considered
  WARNING      alternative in the search space but M* inside the mass range (precursor mass cannot
               separate them for heavier peptides; fragment evidence or mass-error distributions
               must); OR alternative absent but M* above the whole mass range (identity rests on
               precursor tolerance and mass calibration alone)
  FAIL         alternative absent AND heavier peptides fall at or above M*: when observed positive
               peptide masses are given, FAIL needs >= SEARCHSPACE_FAIL_SHARE of them at or above M*;
               otherwise the declared/default mass range is used
"""
from __future__ import annotations

import numpy as np

from ..constants import (DEFAULT_PEPTIDE_MASS_RANGE, FAIL, MASS_H, MASS_O, MASS_P, MASS_S, PASS,
                         SEARCHSPACE_FAIL_SHARE, UNDECIDABLE, WARNING)
from .common import rnd, undecidable

TEST = "search_space"

DELTAS = {"persulfidation": MASS_S, "dioxidation": 2 * MASS_O, "trioxidation": 3 * MASS_O,
          "oxidised_persulfide": MASS_S + MASS_O,
          "phosphorylation": MASS_H + MASS_P + 3 * MASS_O, "sulfation": MASS_S + 3 * MASS_O}
ALIASES = {"persulfidation": "persulfidation", "persulfide": "persulfidation", "sulfhydration": "persulfidation",
           "sulfide": "persulfidation", "s-sulfhydration": "persulfidation",
           "dioxidation": "dioxidation", "sulfinylation": "dioxidation", "sulfinic": "dioxidation",
           "sulfinic_acid": "dioxidation",
           "trioxidation": "trioxidation", "sulfonylation": "trioxidation", "sulfonic": "trioxidation",
           "sulfonic_acid": "trioxidation",
           "oxidised_persulfide": "oxidised_persulfide", "oxidized_persulfide": "oxidised_persulfide",
           "phosphorylation": "phosphorylation", "phospho": "phosphorylation", "phosphoryl": "phosphorylation",
           "sulfation": "sulfation", "sulphation": "sulfation", "sulfo": "sulfation", "sulpho": "sulfation",
           "sulfotyrosine": "sulfation"}
ALTERNATIVES = {"persulfidation": ["dioxidation"], "dioxidation": ["persulfidation"],
                "trioxidation": ["oxidised_persulfide"], "oxidised_persulfide": ["trioxidation"],
                "phosphorylation": ["sulfation"], "sulfation": ["phosphorylation"]}


def canonical(name):
    """Map a modification name to the table key; a trailing residue spec such as '(C)' is ignored."""
    if name is None:
        return None
    base = name.split("(")[0].strip().lower().replace(" ", "_").replace("-", "_")
    return ALIASES.get(base, None)


def resolving_mass(delta_da, ppm):
    return delta_da / (ppm * 1e-6)


def run(ds, cfg):
    mod = cfg.get("modification")
    if cfg.get("identity_readout") == "chemistry_inferred":
        return undecidable(TEST, "identity readout is chemistry_inferred: the modification identity was assigned "
                                 "by the enrichment chemistry, not by a measured mass; no search space confirms it",
                           details={"modification": mod})
    cm = canonical(mod)
    if cm is None or cm not in ALTERNATIVES:
        return undecidable(TEST, f"modification {mod!r} is not covered by the built-in isobaric table",
                           details={"covered": sorted(ALTERNATIVES)})
    ppm = cfg.get("precursor_ppm")
    if ppm is None:
        return undecidable(TEST, "precursor tolerance (--precursor-ppm) not supplied")
    searched = {canonical(m) for m in cfg.get("search_mods", []) if canonical(m)}
    if cm not in searched:
        searched_note = "claimed modification not listed in --search-mods; assumed searched"
    else:
        searched_note = ""
    lo_m, hi_m = cfg.get("mass_range") or DEFAULT_PEPTIDE_MASS_RANGE
    masses = None
    if ds is not None and ds.peptide_mass is not None:
        m = ds.peptide_mass[(ds.label == 1) & np.isfinite(ds.peptide_mass)]
        if m.size:
            masses = m
    rows = []
    worst = PASS
    order = {PASS: 0, WARNING: 1, FAIL: 2}
    for alt in ALTERNATIVES[cm]:
        delta = abs(DELTAS[cm] - DELTAS[alt])
        mstar = resolving_mass(delta, ppm)
        present = alt in searched
        share_above = None
        if masses is not None:
            share_above = float((masses >= mstar).mean())
            heavy = share_above >= SEARCHSPACE_FAIL_SHARE
            inside_range = share_above > 0
        else:
            heavy = inside_range = hi_m >= mstar
        if present:
            st = WARNING if inside_range else PASS
        else:
            st = FAIL if heavy else WARNING
        rows.append({"alternative": alt, "delta_da": rnd(delta, 6), "resolving_mass_da": rnd(mstar, 1),
                     "alternative_in_search_space": present, "share_positive_peptides_at_or_above_resolving_mass":
                     rnd(share_above), "status": st})
        if order[st] > order[worst]:
            worst = st
    head = max(rows, key=lambda r: order[r["status"]])
    reason = (f"{cm} vs {head['alternative']}: delta {head['delta_da']} Da, resolving mass {head['resolving_mass_da']} Da "
              f"at {ppm} ppm; alternative {'in' if head['alternative_in_search_space'] else 'NOT in'} the search space")
    return {"test": TEST, "status": worst, "reason": reason, "statistic": "resolving mass (Da) of the closest isobaric alternative",
            "null": None, "margin": None, "estimate": head["resolving_mass_da"], "ci": None, "p_value": None,
            "n_positive": None if masses is None else int(masses.size), "n_background": None,
            "details": {"modification": cm, "precursor_ppm": ppm, "search_mods_canonical": sorted(searched),
                        "mass_range_used": None if masses is not None else [lo_m, hi_m],
                        "used_observed_peptide_masses": masses is not None, "alternatives": rows,
                        "note": searched_note}}
