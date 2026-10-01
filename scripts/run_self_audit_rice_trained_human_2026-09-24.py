"""Reverse-direction self-audit: the ranker trained on rice PXD072089, applied zero-shot to human PXD044043.

Pre-registered in protocols/self_audit_reverse_rice_trained_preregistration_2026-09-24.json (G0 checks its SHA-256
against results/self_audit_reverse_rice_trained_registered_2026-09-24.sha256). Everything is the human-trained /
rice-audited script (header below) with the cohorts swapped, plus the A2 permutation-distribution null control run
inside this script as the gating null (1,000 draws per cell, seed 20260924), as the pre-registration fixes. The
single-draw NC1a is still computed and reported for continuity, but is not gating here.
Earlier header follows.
Self-audit of the ranker trained on human PXD044043, applied zero-shot to rice PXD072089.

Replaces scripts/run_self_audit_public_cohorts_2026-09-20.py now that the ranker is trained on public data
(scripts/train_v2_public_human_2026-09-23.py) instead of the removed tomato cohort. Differences from that script,
everything else identical:
  - ONE cohort, not two: human_PXD044043 is training data now, not a held-out audit cohort, so it is not scored.
    Rice_PXD072089 is the only genuinely external cohort left.
  - SCORES points at results/v2_public_human_rice_scores_2026-09-23.csv (columns v2h_chem_score, v2h_full_score,
    observed_in_rice), not external/v2_retrospective_site_scores.csv (the removed tomato-trained model's scores).
  - PC1 no longer reproduces stored TOMATO outputs (that data and its derived scores are gone). It instead checks
    that the imported instrument files are BYTE-IDENTICAL to the copies whose correctness was already established
    by PC1 of the 2026-09-20 run (results/self_audit_public_cohorts_2026-09-20_audit.json), i.e. it reuses the
    already-validated instrument rather than re-deriving its correctness from data that no longer exists. This is
    a weaker check than re-deriving correctness (it cannot catch a bug that was already present on 2026-09-20 and
    still is), but the instrument has not been touched since, and re-deriving correctness would need exactly the
    tomato-trained scores this round exists to stop using.
  - NC2 (protein overlap with tomato) is dropped: nothing here can overlap with data that is no longer read.
  - The verdict rule (E1 replicates / E2 partial / E3 clears / E4, from the pre-registration) is UNCHANGED, but with
    one cohort instead of two "both cohorts" collapses to "the one cohort", stated explicitly in the audit output
    rather than silently narrowed.

The pre-registration (scripts/preregister_self_audit_public_2026-09-20.py, protocols/
self_audit_public_cohorts_preregistration_2026-09-20.json) is UNCHANGED and still gates this run: it was written
before ANY public score was read, tomato or human-trained, and nothing in its branches (E1-E4) or its amendment A1
(NC1 split into global/within-protein) refers to which cohort trained the model. Re-using it here is therefore
applying a criterion fixed in advance to a new, not-yet-read dataset, not relaxing or re-writing a criterion after
seeing a result.

INTERPRETER: project venv (Python 3.9.6 / numpy 2.0.2), matching the instrument's own declared interpreter.
"""
import csv, gzip, hashlib, importlib.util, json, pathlib, re, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
PREREG = ROOT / "scripts/preregister_self_audit_public_2026-09-20.py"
PREREG_AUDIT = ROOT / "results/self_audit_public_prereg_2026-09-20_audit.json"
PROTOCOL = ROOT / "protocols/self_audit_public_cohorts_preregistration_2026-09-20.json"
INSTR = ROOT / "scripts/run_self_audit_five_proteases.py"
FLANKP = ROOT / "scripts/run_self_audit_cleavage_flank.py"
PRIOR_AUDIT = ROOT / "results/self_audit_public_cohorts_2026-09-20_audit.json"
SCORES = ROOT / "results/v2_public_rice_human_scores_2026-09-24.csv"
REV_PREREG = ROOT / "protocols/self_audit_reverse_rice_trained_preregistration_2026-09-24.json"
REV_PREREG_SHA = ROOT / "results/self_audit_reverse_rice_trained_registered_2026-09-24.sha256"
OUT_A2 = ROOT / "results/self_audit_rice_trained_human_A2_2026-09-24.csv"
N_PERM, A2_SEED = 1000, 20260924
OUT_CSV = ROOT / "results/self_audit_rice_trained_human_2026-09-24.csv"
OUT_5P = ROOT / "results/self_audit_rice_trained_human_five_proteases_2026-09-24.csv"
OUT_AUDIT = ROOT / "results/self_audit_rice_trained_human_2026-09-24_audit.json"
SELF = pathlib.Path(__file__).read_bytes()

RICE_FASTA = ROOT / "external/proteomes/hsa.fasta.gz"  # the AUDITED cohort's proteome (human); name kept
MODELS = {"v2r_chem_trained_on_rice": "v2r_chem_score", "v2r_full_trained_on_rice": "v2r_full_score"}
PERM_SEED = 20260920


def sha(b):
    return hashlib.sha256(b).hexdigest()


def read_fasta_gz(path):
    out, acc, buf = {}, None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    out[acc] = "".join(buf)
                m = re.match(r">\w+\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf)
    return out


def main():
    if OUT_AUDIT.exists():
        sys.exit("REFUSE: audit exists")
    if sha(REV_PREREG.read_bytes()) != REV_PREREG_SHA.read_text().split()[0]:
        sys.exit("REFUSE: G0 reverse pre-registration changed")
    gates, fails, notes = {}, [], []

    # ---- A0: pre-registration unchanged ----
    pa = json.loads(PREREG_AUDIT.read_text(encoding="utf-8"))
    gates["A0_prereg_recorded"] = pa["script_sha256"]
    gates["A0_prereg_now"] = sha(PREREG.read_bytes())
    gates["A0_protocol_recorded"] = pa["protocol_sha256"]
    gates["A0_protocol_now"] = sha(PROTOCOL.read_bytes())
    gates["A0_unchanged"] = (gates["A0_prereg_recorded"] == gates["A0_prereg_now"]
                             and gates["A0_protocol_recorded"] == gates["A0_protocol_now"])
    if not gates["A0_unchanged"]:
        sys.exit("REFUSED: the pre-registration changed after it was written")

    # ---- PC1 (amended): the instrument files are byte-identical to the copies validated on 2026-09-20 ----
    prior = json.loads(PRIOR_AUDIT.read_text(encoding="utf-8"))
    now_instr, now_flank = sha(INSTR.read_bytes()), sha(FLANKP.read_bytes())
    gates["PC1_instrument_unchanged_since_validated_run"] = dict(
        instrument=dict(recorded=prior["gates"]["instrument"]["sha256"], now=now_instr,
                        matches=now_instr == prior["gates"]["instrument"]["sha256"]),
        flank=dict(recorded=prior["gates"]["instrument_flank"]["sha256"], now=now_flank,
                   matches=now_flank == prior["gates"]["instrument_flank"]["sha256"]),
        prior_run_PC1_pass=prior["gates"]["PC1_pass"])
    gates["PC1_pass"] = (gates["PC1_instrument_unchanged_since_validated_run"]["instrument"]["matches"]
                         and gates["PC1_instrument_unchanged_since_validated_run"]["flank"]["matches"]
                         and prior["gates"]["PC1_pass"])
    if not gates["PC1_pass"]:
        OUT_AUDIT.write_text(json.dumps(dict(gates=gates, gates_failed=["PC1"]), ensure_ascii=False, indent=2))
        sys.exit("REFUSED: PC1 -- instrument changed since its correctness was last validated")

    # ---- import the instrument (unchanged code, confirmed above) ----
    assert "__main__" in INSTR.read_text(encoding="utf-8"), "instrument not import-safe"
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("fivep", INSTR)
    F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    spec2 = importlib.util.spec_from_file_location("flankmod", FLANKP)
    FL = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(FL)
    gates["instrument"] = dict(name=INSTR.name, sha256=now_instr, DISTAL=list(F.DISTAL), PROXIMAL=F.PROXIMAL,
                               REPLICATES=F.REPLICATES, SEED=F.SEED, proteases=list(F.PROTEASES))
    gates["instrument_flank"] = dict(name=FLANKP.name, sha256=now_flank,
                                     SFE006_KR_OFFSETS=list(FL.SFE006_KR_OFFSETS))

    def interval(attr, score, groups):
        fn = F.make_weighted_auc(attr, score)
        point = fn(np.ones(len(attr), dtype=float))
        lo_hi, _ = F.bootstrap_interval(fn, groups)
        return float(point), float(lo_hi[0]), float(lo_hi[1])

    def attrs_for(seqs, accs, poss, rule):
        d, p = [], []
        for a, q in zip(accs, poss):
            x, y, _ = F.band_flags(seqs[a], q, rule)
            d.append(x); p.append(y)
        return np.asarray(d, int), np.asarray(p, int)

    tryp = F.PROTEASES["Trypsin"]

    # ---- rice cohort ----
    seqs = read_fasta_gz(RICE_FASTA)
    allrows = list(csv.DictReader(open(SCORES, encoding="utf-8-sig")))
    kept = [r for r in allrows
            if r["protein"] in seqs
            and 0 < int(r["position"]) <= len(seqs[r["protein"]])
            and seqs[r["protein"]][int(r["position"]) - 1] == "C"]
    accs = [r["protein"] for r in kept]
    poss = [int(r["position"]) for r in kept]
    y = np.asarray([int(float(r["observed_in_human"])) for r in kept])
    g = np.asarray(accs)
    cov = dict(rows_in_table=len(allrows), rows_kept=len(kept),
              proteins_in_table=len({r["protein"] for r in allrows}),
              proteins_matched=len(set(accs)), positives=int(y.sum()))
    gates["coverage_human_PXD044043"] = cov
    if cov["rows_kept"] < cov["rows_in_table"]:
        notes.append(f"human_PXD044043: {cov['rows_in_table'] - cov['rows_kept']} of {cov['rows_in_table']} sites "
                     f"dropped because the protein is absent from the reference proteome or the recorded position "
                     f"is not a cysteine")

    d, p = attrs_for(seqs, accs, poss, tryp)
    fk = np.asarray([FL.flank_flag(seqs[a], q, FL.SFE006_KR_OFFSETS) or 0 for a, q in zip(accs, poss)], dtype=int)
    gates["PC3_human_distal_rate"] = round(float(d.mean()), 4)
    gates["PC3_human_proximal_rate"] = round(float(p.mean()), 4)
    if not (0.10 < d.mean() < 0.99):
        fails.append(f"PC3 judged negative: human distal rate {d.mean():.4f} outside (0.10, 0.99)")

    lpt, llo, lhi = interval(d, y.astype(float), g)
    gates["PC2_human_label_distal_auc"] = dict(auc=round(lpt, 4), ci=[round(llo, 4), round(lhi, 4)])

    rows_out, rows_5p, cohort_res = [], [], {}
    for mname, col in MODELS.items():
        s = np.asarray([float(r[col]) for r in kept])
        for band, attr in (("distal_6_12", d), ("proximal_le_5", p), ("kr_sfe006_offsets", fk)):
            for stratum, mask in (("all", np.ones(len(y), bool)), ("label_positive", y == 1),
                                  ("label_negative", y == 0)):
                if mask.sum() < 20 or len(np.unique(attr[mask])) < 2:
                    continue
                pt, lo, hi = interval(attr[mask], s[mask], g[mask])
                rows_out.append(dict(cohort="human_PXD044043", model=mname, band=band, stratum=stratum,
                                     n=int(mask.sum()), attribute_rate=round(float(attr[mask].mean()), 4),
                                     auc=round(pt, 4), ci_low=round(lo, 4), ci_high=round(hi, 4)))
                cohort_res[(mname, band, stratum)] = (pt, lo, hi)

        # NC1a (global permutation, a true null, GATING) / NC1b (within-protein, reported, not gating) -- as
        # amended 2026-09-20 (A1), unchanged here
        rng = np.random.default_rng(PERM_SEED)
        sp = s.copy()
        order = np.argsort(g, kind="stable")
        gs = g[order]
        starts = np.flatnonzero(np.concatenate([[True], gs[1:] != gs[:-1]]))
        for b in np.split(order, starts[1:]):
            sp[b] = s[rng.permutation(b)]
        rng_g = np.random.default_rng(PERM_SEED + 1)
        for stratum, mask in (("label_positive", y == 1), ("label_negative", y == 0)):
            sg = s[mask].copy()
            rng_g.shuffle(sg)
            pt, lo, hi = interval(d[mask], sg, g[mask])
            covers = bool(lo <= 0.5 <= hi)
            gates[f"NC1a_global_human_{mname}_{stratum}"] = dict(auc=round(pt, 4), ci=[round(lo, 4), round(hi, 4)],
                                                                covers_half=covers)
            if False and not covers:  # not gating here: A2 below is the pre-registered gating null
                fails.append(f"NC1a judged negative: global permutation still separates the attribute in "
                             f"human/{mname}/{stratum}: {pt:.4f} [{lo:.4f}, {hi:.4f}]. Per the protocol the "
                             f"estimator is broken and no result from this round may be used.")
            pt2, lo2, hi2 = interval(d[mask], sp[mask], g[mask])
            gates[f"NC1b_within_protein_human_{mname}_{stratum}"] = dict(
                auc=round(pt2, 4), ci=[round(lo2, 4), round(hi2, 4)], covers_half=bool(lo2 <= 0.5 <= hi2),
                reads_as="between-protein floor: the part of the pooled association that survives destroying "
                        "all within-protein ordering")
    per_cohort = {"human_PXD044043": cohort_res}

    for pname, rule in F.PROTEASES.items():
        dd, pp = attrs_for(seqs, accs, poss, rule)
        for mname, col in MODELS.items():
            s = np.asarray([float(r[col]) for r in kept])
            for band, attr in (("distal_6_12", dd), ("proximal_le_5", pp)):
                for stratum, mask in (("all", np.ones(len(y), bool)), ("label_positive", y == 1),
                                      ("label_negative", y == 0)):
                    if mask.sum() < 20 or len(np.unique(attr[mask])) < 2:
                        continue
                    pt, lo, hi = interval(attr[mask], s[mask], g[mask])
                    rows_5p.append(dict(cohort="human_PXD044043", protease=pname, model=mname, band=band,
                                        stratum=stratum, n=int(mask.sum()),
                                        attribute_rate=round(float(attr[mask].mean()), 4),
                                        auc=round(pt, 4), ci_low=round(lo, 4), ci_high=round(hi, 4)))

    # ---- A2 (pre-registered, gating): 1,000-draw permutation distribution per cell ----
    rng_a2 = np.random.default_rng(A2_SEED)
    a2_rows = []
    for mname, col in MODELS.items():
        s_all = np.asarray([float(r[col]) for r in kept])
        for stratum in ("label_positive", "label_negative"):
            mask = (y == 1) if stratum == "label_positive" else (y == 0)
            a_, sc = d[mask], s_all[mask]
            ones = np.ones(len(a_))
            observed = F.make_weighted_auc(a_, sc)(ones)
            null = np.asarray([F.make_weighted_auc(a_, rng_a2.permutation(sc))(ones) for _ in range(N_PERM)])
            lo_, hi_ = np.percentile(null, [2.5, 97.5])
            covers = bool(lo_ <= 0.5 <= hi_)
            if not covers:
                fails.append(f"A2 null fails in {mname}/{stratum}: [{lo_:.4f}, {hi_:.4f}]")
            a2_rows.append(dict(cohort="human_PXD044043", model=mname, stratum=stratum, n=int(mask.sum()),
                                observed_auc=round(float(observed), 4), null_mean=round(float(null.mean()), 4),
                                null_p2_5=round(float(lo_), 4), null_p97_5=round(float(hi_), 4),
                                null_covers_half=covers,
                                permutation_p=round((1 + int((null >= observed).sum())) / (N_PERM + 1), 4),
                                n_permutations=N_PERM))
    with open(OUT_A2, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(a2_rows[0])); w.writeheader(); w.writerows(a2_rows)
    gates["A2_all_cells_cover_half"] = all(r["null_covers_half"] for r in a2_rows)

    # ---- verdict, pre-declared rule, one cohort ----
    above = {}
    for stratum in ("label_positive", "label_negative"):
        k = ("v2r_chem_trained_on_rice", "distal_6_12", stratum)
        if k in cohort_res:
            above[stratum] = cohort_res[k][1] > 0.5
    all_above = all(above.values()) and len(above) == 2
    none_above = not any(above.values())
    if all_above:
        verdict, why = "E1_replicates", "interval above 0.5 in both label strata (human audited; rice trains the model)"
    elif none_above:
        verdict, why = "E3_clears", "no within-stratum interval lies above 0.5: the rice-trained model does " \
                                    "not rank human cysteines by cleavage geometry out-of-sample"
    else:
        verdict, why = "E2_partial", "strata disagree"
    gates["verdict"] = verdict
    gates["verdict_why"] = why
    gates["verdict_inputs"] = {s: bool(v) for s, v in above.items()}
    gates["single_cohort_declared"] = "reverse direction: rice trains, human is audited"

    for path, rows in ((OUT_CSV, rows_out), (OUT_5P, rows_5p)):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    train_audit = json.loads((ROOT / "results/train_v2_public_rice_2026-09-24_audit.json").read_text())
    OUT_AUDIT.write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=sha(SELF), interpreter=sys.version.split()[0],
        numpy=np.__version__, training_run=dict(script="scripts/train_v2_public_rice_2026-09-24.py",
                                                script_sha256=train_audit["script_sha256"],
                                                train_report=train_audit["report"]),
        preregistration=dict(script=PREREG.name, protocol=PROTOCOL.name,
                             sha256_verified_unchanged=gates["A0_unchanged"]),
        declared_deviation="bootstrap clusters by PROTEIN, not homology component: the rice_PXD072089 training "
                           "cohort has no homology clustering (components == proteins; scripts/"
                           "build_public_refit_inputs_2026-09-21.py's own note). Cross-validation folds during "
                           "training were therefore protein-disjoint but not verified homology-disjoint, and the "
                           "bootstrap here clusters by protein for the same reason. Intervals are narrower than "
                           "true homology-component clustering would give and must be quoted with this note.",
        gates=gates, gates_failed=fails, notes=notes,
        what_this_does_NOT_establish=[
            "whether the model is a good persulfidation predictor",
            "anything about tomato, or about any chemistry other than persulfidation",
            "whether cross-species transfer (rice-trained, human-tested) differs from same-species transfer, "
            "which was not tested",
        ]), ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(dict(verdict=verdict, why=why, pc1=gates["PC1_pass"], rows=len(rows_out),
                          rows_5p=len(rows_5p), gates_failed=fails), ensure_ascii=False))


if __name__ == "__main__":
    main()
