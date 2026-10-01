"""Read the verdicts for the public refit of manuscript 3.3 and 3.4.

Criteria are NOT defined here. They were fixed in scripts/preregister_public_refit_2026-09-21.py and
protocols/public_refit_3_3_and_3_4_preregistration_2026-09-21.json, written before any arm was
fitted. Gate A0 refuses to run if either changed. Branch P2 -- the outcome that CONTRADICTS the
manuscript -- was fixed there in advance.

The refit itself ran on amax and returned out-of-fold SCORES ONLY. Every verdict statistic is
computed HERE, by the already-registered instrument, imported not reimplemented (CLAUDE.md 9.24.2):
  scripts/run_self_audit_five_proteases.py   DISTAL, band_flags, TOP_K, make_weighted_table,
                                             log_odds_ratio, bootstrap_interval (via its imports)
  scripts/common.py                          protein_rows (the per-protein AUC that defines
                                             `within_protein_auc_mean`)
so the numbers are produced by the same code that produced the tomato ablation's 3.7116.

INTERPRETER: project venv (Python 3.9.6 / numpy 2.0.2).
"""
import csv, gzip, hashlib, importlib.util, json, pathlib, re, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
PREREG = ROOT / "scripts/preregister_public_refit_2026-09-21.py"
PREREG_AUDIT = ROOT / "results/public_refit_prereg_2026-09-21_audit.json"
PROTOCOL = ROOT / "protocols/public_refit_3_3_and_3_4_preregistration_2026-09-21.json"
OUTDIR = ROOT / "external/public_refit_out_2026-09-21"
INSTR = ROOT / "scripts/run_self_audit_five_proteases.py"
OUT_CSV = ROOT / "results/public_refit_ablation_2026-09-21.csv"
OUT_DET = ROOT / "results/public_refit_detectability_only_2026-09-21.csv"
OUT_AUDIT = ROOT / "results/public_refit_ablation_2026-09-21_audit.json"
SELF = pathlib.Path(__file__).read_bytes()

COHORTS = {"human_PXD044043": ROOT / "external/proteomes/hsa.fasta.gz",
           "rice_PXD072089": ROOT / "external/proteomes/osa.fasta.gz"}
PRIMARY = "human_PXD044043"
ARMS_33 = ("control", "drop_w710", "drop_all_basic")
BOOT, SEED = 5000, 20260921


def sha(b):
    return hashlib.sha256(b).hexdigest()


def read_oof(cohort, arm):
    f = OUTDIR / cohort / f"{cohort}_{arm}_oof.csv"
    rows = list(csv.DictReader(open(f, encoding="utf-8-sig")))
    return dict(accession=np.asarray([r["accession"] for r in rows]),
                position=np.asarray([int(r["position"]) for r in rows]),
                y=np.asarray([int(r["label"]) for r in rows]),
                score=np.asarray([float(r["blend"]) for r in rows]),
                file_sha256=sha(f.read_bytes()), n=len(rows))


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
    gates, fails, notes = {}, [], []

    # ---- A0 -----------------------------------------------------------------
    pa = json.loads(PREREG_AUDIT.read_text(encoding="utf-8"))
    gates["A0_unchanged"] = (pa["script_sha256"] == sha(PREREG.read_bytes())
                             and pa["protocol_sha256"] == sha(PROTOCOL.read_bytes()))
    gates["A0_prereg_sha256"] = pa["script_sha256"]
    if not gates["A0_unchanged"]:
        print("REFUSED: the pre-registration changed after it was written"); sys.exit(1)

    # ---- import the registered instrument ------------------------------------
    assert "__main__" in INSTR.read_text(encoding="utf-8")
    spec = importlib.util.spec_from_file_location("fivep", INSTR)
    F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    from common import protein_rows
    gates["instrument"] = dict(name=INSTR.name, sha256=sha(INSTR.read_bytes()),
                               DISTAL=list(F.DISTAL), TOP_K=F.TOP_K,
                               REPLICATES=F.REPLICATES, SEED=F.SEED)

    tryp = F.PROTEASES["Trypsin"]
    rows_out, rows_det, verdicts = [], [], {}
    for cohort, fa in COHORTS.items():
        seqs = read_fasta_gz(fa)
        base = read_oof(cohort, "control")
        acc, pos = base["accession"], base["position"]
        groups = acc
        keep = np.asarray([a in seqs and 0 < q <= len(seqs[a]) for a, q in zip(acc, pos)])
        gates[f"{cohort}_rows_with_sequence"] = dict(kept=int(keep.sum()), total=int(len(keep)))
        dist = np.asarray([F.band_flags(seqs[a], q, tryp)[0] if k else -1
                           for a, q, k in zip(acc, pos, keep)], dtype=int)

        per_arm = {}
        for arm in ARMS_33 + ("label_permuted", "detect_only", "full"):
            d = read_oof(cohort, arm)
            if not (np.array_equal(d["accession"], acc) and np.array_equal(d["position"], pos)):
                fails.append(f"row order differs for {cohort}/{arm}")
                continue
            y_arm, s, m = d["y"], d["score"], keep
            order = np.argsort(-s[m], kind="stable")
            in_top = np.zeros(int(m.sum()), dtype=int)
            in_top[order[:F.TOP_K]] = 1
            table = F.make_weighted_table(in_top, dist[m])
            lor = F.log_odds_ratio(*table(np.ones(int(m.sum()))))
            liv, _ = F.bootstrap_interval(lambda w: F.log_odds_ratio(*table(w)), groups[m])
            pr = protein_rows(y_arm[m], s[m], acc[m], acc[m])
            wp = float(np.mean([r["auc"] for r in pr]))
            per_arm[arm] = dict(lor=float(lor), lor_ci=[float(liv[0]), float(liv[1])],
                                rate_top=float(dist[m][order[:F.TOP_K]].mean()),
                                rate_rest=float(dist[m][order[F.TOP_K:]].mean()),
                                wp_auc=wp, per_protein={r["accession"]: r["auc"] for r in pr})
            rows_out.append(dict(cohort=cohort, arm=arm, n=int(m.sum()),
                                 top100_distal_rate=round(per_arm[arm]["rate_top"], 4),
                                 rate_rest=round(per_arm[arm]["rate_rest"], 4),
                                 log2_or=round(lor, 4), ci_low=round(liv[0], 4), ci_high=round(liv[1], 4),
                                 within_protein_auc=round(wp, 4), n_proteins_scored=len(pr)))
        if fails:
            continue

        # ---- NC1: a permuted label must give neither ranking nor enrichment
        nc1 = per_arm["label_permuted"]
        crosses = nc1["lor_ci"][0] <= 0.0 <= nc1["lor_ci"][1]
        rng = np.random.default_rng(SEED)
        vals = np.asarray(list(nc1["per_protein"].values()))
        bs = np.asarray([vals[rng.integers(0, len(vals), len(vals))].mean() for _ in range(BOOT)])
        wp_ci = (float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)))
        covers = wp_ci[0] <= 0.5 <= wp_ci[1]
        gates[f"NC1_{cohort}"] = dict(top100_lor=round(nc1["lor"], 4),
                                      top100_lor_ci=[round(x, 4) for x in nc1["lor_ci"]],
                                      crosses_zero=bool(crosses), wp_auc=round(nc1["wp_auc"], 4),
                                      wp_auc_ci=[round(wp_ci[0], 4), round(wp_ci[1], 4)],
                                      wp_covers_half=bool(covers))
        if not (crosses and covers):
            fails.append(f"NC1 judged negative for {cohort}: {gates[f'NC1_{cohort}']}. "
                         f"Per the protocol nothing from this round is usable.")
        gates[f"PC2_{cohort}_control_wp_auc"] = round(per_arm["control"]["wp_auc"], 4)
        if cohort == PRIMARY and per_arm["control"]["wp_auc"] <= 0.5:
            fails.append("PC2 judged negative: control within-protein AUC <= 0.5 in the primary cohort")

        # ---- 3.3 branch, by the pre-declared rule
        c = per_arm["control"]["lor_ci"]
        ablated = {a: per_arm[a]["lor_ci"] for a in ("drop_w710", "drop_all_basic")}
        if c[0] <= 0.0 <= c[1]:
            branch = "P3_no_baseline"
        elif all(not (v[0] <= 0.0 <= v[1]) for v in ablated.values()):
            branch = "P1_survives"
        else:
            branch = "P2_removed"
        verdicts[cohort] = dict(branch_3_3=branch, control=c, ablated=ablated)

        # ---- 3.4: recovery ratio, paired protein bootstrap
        pf, pt = per_arm["full"]["per_protein"], per_arm["detect_only"]["per_protein"]
        cp = sorted(set(pf) & set(pt))
        vf = np.asarray([pf[p] for p in cp]); vt = np.asarray([pt[p] for p in cp])
        rng = np.random.default_rng(SEED)

        def ratio(a, b):
            den = a.mean() - 0.5
            return (b.mean() - 0.5) / den if den > 0 else np.nan

        point = ratio(vf, vt)
        draws = []
        for _ in range(BOOT):
            i = rng.integers(0, len(cp), len(cp))
            r = ratio(vf[i], vt[i])
            if np.isfinite(r):
                draws.append(r)
        lo, hi = float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))
        rows_det.append(dict(cohort=cohort, n_proteins=len(cp),
                             wp_auc_full=round(float(vf.mean()), 4),
                             wp_auc_detect_only=round(float(vt.mean()), 4),
                             recovery_ratio=round(float(point), 4),
                             ci_low=round(lo, 4), ci_high=round(hi, 4),
                             branch_3_4="Q1" if lo > 0.50 else "Q2"))
        verdicts[cohort]["branch_3_4"] = rows_det[-1]["branch_3_4"]

    if not rows_out or not rows_det:
        print("REFUSED:", fails); sys.exit(1)
    for path, rows in ((OUT_CSV, rows_out), (OUT_DET, rows_det)):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    gates["verdicts"] = verdicts
    if PRIMARY in verdicts:
        gates["PRIMARY_verdict_3_3"] = verdicts[PRIMARY]["branch_3_3"]
        gates["PRIMARY_verdict_3_4"] = verdicts[PRIMARY]["branch_3_4"]
    OUT_AUDIT.write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=sha(SELF),
        interpreter=sys.version.split()[0], numpy=np.__version__,
        preregistration=dict(script=PREREG.name, protocol=PROTOCOL.name, sha256_verified_unchanged=True),
        bootstrap=dict(replicates=BOOT, seed=SEED, unit="protein (declared deviation from house rule 2)"),
        gates=gates, gates_failed=fails, notes=notes,
        what_this_does_NOT_establish=[
            "it does not restore 3.3 or 3.4 as statements about the deployed tomato model",
            "it says nothing about whether the model is a good persulfidation predictor",
            "both cohorts are persulfidation; nothing transfers to other chemistries",
            "a null result on rice is not evidence of absence (817 positives, 640 proteins)",
            "no structural features for external sequences: absolute performance is not comparable "
            "across cohorts, only arm-to-arm within a cohort",
        ]), ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(dict(primary_3_3=gates.get("PRIMARY_verdict_3_3"),
                          primary_3_4=gates.get("PRIMARY_verdict_3_4"),
                          per_cohort={c: (v["branch_3_3"], v["branch_3_4"]) for c, v in verdicts.items()},
                          gates_failed=fails), ensure_ascii=False))


if __name__ == "__main__":
    main()
