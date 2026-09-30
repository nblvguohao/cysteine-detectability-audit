"""Figures 2-5 of manuscript R22, redrawn in Nature Communications house style. Figure 1 is the new overview
(scripts/plot_fig1_overview_2026-09-21.py).

RENUMBERING (the journal: "Figures should be numbered ... in the order of occurrence in the text"; in R21 the old F4 was
cited first):  Fig. 2 = old F4 (search space) | Fig. 3 = old F1 (three axes) | Fig. 4 = old F2 (claim re-tests) |
Fig. 5 = old F3 (self-audit on public cohorts).

WHAT CHANGED AND WHY (user: composite panels squashed the axes; learn from published NC figures)
  * drawn at print size, 183 mm wide, Arial 7 pt, bold 8 pt panel letters, every line >= 1 pt (journal text);
  * each panel gets the height its row count needs; no panel is shorter than 45 mm;
  * explanations that were printed inside the plots (subtitles, dataset notes, verdict words, row tags, the
    no-baseline list) move to the figure legends, as published NC figures do; plots keep only labels and values;
  * one colour per verdict, the same in Figs 3 and 4; Okabe-Ito colours for cohorts.

DATA: NOTHING IS RE-READ BY NEW CODE WHERE A READER ALREADY EXISTS. Figs 3 and 4 call the collectors of the original plot
scripts (plot_manuscript_f1_three_axes.collect_panel_a/b/c and plot_manuscript_f2_claim_retest.collect), whose sha256
is recorded; Fig. 2 re-applies the original script's xcorr binning (copied verbatim below, declared); Fig. 5 reads the
same three public tables as scripts/plot_manuscript_f3_public_2026-09-21.py. Every value drawn is logged with its raw
stored string and written to a Source Data CSV per figure.

GATES (all four figures): G1 every drawn value equals float(raw stored cell); G2 minimum line width >= 1.0 pt;
G3 only font sizes 7 and 8; G4 zero text collisions (tree detector, positive control); G5 no drawn text outside the
canvas; G6 no panel axes shorter than 45 mm.
INTERPRETER: /Users/lyuguohao/Documents/lnrna复现/.venv/bin/python
"""
import csv, hashlib, importlib.util, json, pathlib, sys
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures_r22"
SD = ROOT / "source_data_r22"
SELF = pathlib.Path(__file__).read_bytes()


def load(name):
    p = ROOT / "scripts" / name
    assert name.startswith("nc_figure_style") or "__main__" in p.read_text(), name
    spec = importlib.util.spec_from_file_location(name[:-3].replace("-", "_"), p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m, hashlib.sha256(p.read_bytes()).hexdigest()


F1, F1_SHA = load("plot_manuscript_f1_three_axes.py")
F2, F2_SHA = load("plot_manuscript_f2_claim_retest.py")
S, S_SHA = load("nc_figure_style_2026-09-21.py")
CHK, CHK_SHA = load("check_figure_text_overlaps_2026-09-19.py")
S.apply()                       # after the imports: the old modules set their own rcParams on import
# Keep text editable in vector exports when this module is run on its own.
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
width_mm = 183

VC = {"survives": S.BLUE, "attenuated": S.SKY, "null_broken_by_control": S.PINK, "vanishes": S.ORANGE,
      "baseline_contradicts_claim": S.VERMIL, "undecidable": S.GREY}
VT = {"survives": "survives", "attenuated": "attenuated", "null_broken_by_control": "null broken by control",
      "vanishes": "vanishes", "baseline_contradicts_claim": "baseline contradicts claim", "undecidable": "undecidable"}
HUMAN, RICE = S.BLUE, S.ORANGE


def read(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig", newline="")))


class ValueLedger:
    def __init__(self, fig):
        self.fig, self.rows = fig, []

    def add(self, panel, row, field, raw, value, table):
        assert float(raw) == value, (panel, row, field, raw, value)
        self.rows.append(dict(figure=self.fig, panel=panel, row=row, field=field, value=raw, source_table=table))


def lw_ci(ax, y, p, lo, hi, c, marker="o", filled=True, ms=3.6):
    S.ci(ax, y, p, lo, hi, colour=c, marker=marker, filled=filled, ms=ms)


def letter(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=S.FS_PANEL, fontweight="bold", ha="right", va="bottom")


# ------------------------------------------------------------------------------------------------ Fig. 2 (old F4)
def fig2():
    log = ValueLedger("Fig2")
    ties = read(RES / "pxd015307_posthoc_score_ties.csv")
    summary = read(RES / "pxd015307_research_summary.csv")
    both = [r for r in ties if r["chemistries_at_best_xcorr"] == "dioxidation+sulfide"]

    def binx(value):                          # verbatim from scripts/plot_manuscript_f3_f4.py figure_f4()
        value = float(value)
        for threshold in (2.0, 1.5, 1.0, 0.5):
            if value >= threshold:
                return f">= {threshold}"
        return "< 0.5"
    order = ["< 0.5", ">= 0.5", ">= 1.0", ">= 1.5", ">= 2.0"]
    shown = ["<0.5", "0.5-1", "1-1.5", "1.5-2", "≥2"]
    counts = Counter(binx(r["best_xcorr"]) for r in both)

    fig = plt.figure(figsize=(width_mm * S.MM, 62 * S.MM))
    gs = fig.add_gridspec(1, 3, left=0.075, right=0.985, bottom=0.2, top=0.88, wspace=0.55, width_ratios=[1, 1, 1.15])
    ax = fig.add_subplot(gs[0])
    ax.bar(range(5), [counts.get(k, 0) for k in order], color=S.INK, width=0.66, lw=0)
    for i, k in enumerate(order):
        ax.text(i, counts.get(k, 0) + 18, str(counts.get(k, 0)), ha="center", va="bottom")
        log.add("a", k, "n_spectra_tied", str(counts.get(k, 0)), counts.get(k, 0), "pxd015307_posthoc_score_ties.csv (binned)")
    ax.set_xticks(range(5)); ax.set_xticklabels(shown)
    ax.set_ylim(0, max(counts.values()) * 1.15)
    ax.set_xlabel("Best xcorr of the spectrum"); ax.set_ylabel("Spectra with a tie")
    ax.set_title("Two chemistries, one best score"); letter(ax, "a", -0.28)

    ax = fig.add_subplot(gs[1])
    spans = sorted(float(r["ppm_span"]) for r in both)
    ax.plot(spans, [i / len(spans) for i in range(len(spans))], color=S.INK, lw=1.2)
    ax.set_xlabel("Precursor-error span among\ntied interpretations (ppm)"); ax.set_ylabel("Cumulative fraction")
    ax.set_ylim(0, 1.02); ax.set_xlim(0, None)
    ax.set_title("Mass errors far apart"); letter(ax, "b", -0.28)
    log.add("b", "max", "ppm_span", str(max(float(r["ppm_span"]) for r in both)), spans[-1], "pxd015307_posthoc_score_ties.csv")

    ax = fig.add_subplot(gs[2])
    names = [r["arm"].replace("Fig4-L_", "") for r in summary]
    psms = [int(r["n_psms_at_1pct_fdr"]) for r in summary]
    ax.barh(range(len(names)), psms, color=S.GREY, height=0.62, lw=0)
    for i, (n, r) in enumerate(zip(psms, summary)):
        ax.text(n + 6, i, str(n), va="center")
        log.add("c", r["arm"], "n_psms_at_1pct_fdr", r["n_psms_at_1pct_fdr"], n, "pxd015307_research_summary.csv")
        log.add("c", r["arm"], "n_plus32_cys_psms", r["n_plus32_cys_psms"], float(r["n_plus32_cys_psms"]), "pxd015307_research_summary.csv")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names); ax.invert_yaxis()
    ax.set_xlim(0, max(psms) * 1.3); ax.set_xlabel("PSMs at 1% FDR")
    ax.set_title("No +32 identification"); letter(ax, "c", -0.42)
    return fig, log


# ------------------------------------------------------------------------------------------------ Fig. 3 (old F1)
def fig3():
    log = ValueLedger("Fig3")
    pa, pb, pc = F1.collect_panel_a(), F1.collect_panel_b(), F1.collect_panel_c()
    fig = plt.figure(figsize=(width_mm * S.MM, 162 * S.MM))
    gs = fig.add_gridspec(2, 2, left=0.19, right=0.985, bottom=0.1, top=0.96, hspace=0.62, wspace=0.62,
                          width_ratios=[1.25, 1], height_ratios=[1.25, 1])

    # a: share of reachable discriminative power carried by visibility features, per dataset and negative set
    ax = fig.add_subplot(gs[0, 0])
    y, ticks, labels = 0, [], []
    XL = (-1.2, 1.15)
    for rec in pa:
        a, b = rec["A"], rec["B"]
        if a and a["defined"]:
            lw_ci(ax, y + 0.15, a["share"], a["lo"], a["hi"], S.INK, "o", True)
            for f, k in (("share_visibility", "share"), ("lo", "lo"), ("hi", "hi")):
                log.add("a", rec["dataset_id"] + "|NEG_A", f, a[k + "__raw"], a[k], "ptm_detectability_share.csv")
        if b and b["defined"]:
            lo = max(b["lo"], XL[0])
            lw_ci(ax, y - 0.15, b["share"], lo, b["hi"], S.BLUE, "s", False)
            if b["lo"] < XL[0]:
                ax.annotate("", xy=(XL[0], y - 0.15), xytext=(XL[0] + 0.12, y - 0.15),
                            arrowprops=dict(arrowstyle="-|>", color=S.BLUE, lw=S.LW, mutation_scale=6))
            for f, k in (("share_visibility", "share"), ("lo", "lo"), ("hi", "hi")):
                log.add("a", rec["dataset_id"] + "|NEG_B", f, b[k + "__raw"], b[k], "ptm_detectability_share.csv")
        ticks.append(y); labels.append(rec["label"].replace(" (", "\n("))
        y -= 1
    ax.axvline(0, color=S.INK, lw=S.LW); ax.axvline(1, color=S.GREY, lw=S.LW, ls=(0, (2, 2)))
    ax.set_yticks(ticks); ax.set_yticklabels(labels, linespacing=1.25); ax.tick_params(axis="y", length=0)
    ax.set_xlim(*XL); ax.set_ylim(y + 0.4, 0.6)
    ax.set_xlabel("Share of discriminative power carried\nby visibility features alone")
    ax.set_title("Negative-set construction"); letter(ax, "a", -0.02)
    ax.legend(handles=[Line2D([], [], marker="o", color=S.INK, lw=S.LW, ms=3.6, label="not observed (NEG_A)"),
                       Line2D([], [], marker="s", color=S.BLUE, mfc="white", lw=S.LW, ms=3.6, label="observed, unmodified (NEG_B)")],
              loc="upper left", bbox_to_anchor=(-0.25, -0.2), ncol=2, columnspacing=1.0)

    # b: PERS-008 against two backgrounds
    ax = fig.add_subplot(gs[0, 1])
    y, ticks, labels = 0, [], []
    group_names = {"primary": "Whole proteome", "zf_background": "Authors' eligible list"}
    for blk in pb:
        ax.text(-0.02, y + 0.55, "%s (n = %s)" % (group_names[blk["caliber"]], f'{int(blk["n_observations"]):,}'),
                transform=ax.get_yaxis_transform(), fontweight="bold", va="center", ha="right", clip_on=False)
        for est in blk["estimates"]:
            c = {"baseline": S.INK, "after matching": VC[blk["verdict"]], "random control": S.GREY}[est["short"]]
            mk = {"baseline": "o", "after matching": "o", "random control": "x"}[est["short"]]
            lw_ci(ax, y, est["point"], est["lo"], est["hi"], c, mk, est["short"] != "baseline")
            for fld, k in ((est["point_field"], "point"), (est["lo_field"], "lo"), (est["hi_field"], "hi")):
                log.add("b", blk["caliber"], fld, est[k + "__raw"], est[k], "phase2d_claim_retest.csv")
            ticks.append(y); labels.append(est["short"]); y -= 1
        y -= 1.2
    ax.axvline(0, color=S.INK, lw=S.LW)
    ax.set_xlim(-3.2, 6.2); ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 1.6, 1.3)
    ax.set_xlabel("log2 odds ratio (PERS-008)")
    ax.set_title("Background choice"); letter(ax, "b", -0.02)

    # c: five structural claims, attribute rate in positives vs negatives; the precision tier sits in a text column
    ax = fig.add_subplot(gs[1, :])
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, 0.40, pos.height])
    ticks, labels = [], []
    for i, r in enumerate(pc):
        y = -i
        c = VC[r["verdict"]]
        ax.plot([r["neg"], r["pos"]], [y, y], color=c, lw=S.LW, zorder=1)
        ax.plot([r["pos"]], [y], "o", color=c, ms=4.2, zorder=3)
        ax.plot([r["neg"]], [y], "o", mfc="white", mec=c, mew=S.LW, ms=4.2, zorder=3)
        log.add("c", r["claim_id"], "attribute_rate_positive", r["pos__raw"], r["pos"], "phase2e_claim_retest.csv")
        log.add("c", r["claim_id"], "attribute_rate_negative", r["neg__raw"], r["neg"], "phase2e_claim_retest.csv")
        ticks.append(y); labels.append(r["claim_id"])
        ax.text(1.06, y, r["tier"], transform=ax.get_yaxis_transform(), va="center", color=S.INK, linespacing=1.15)
    ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 0.72); ax.set_ylim(-len(pc) + 0.4, 0.6)
    ax.set_xlabel("Fraction of sites carrying the claimed attribute")
    ax.text(1.06, 1.02, "How precisely the attribute was stated", transform=ax.transAxes, fontweight="bold", va="bottom")
    ax.set_title("Attribute definition"); letter(ax, "c", -0.02)
    ax.legend(handles=[Line2D([], [], marker="o", color=S.INK, lw=0, ms=4, label="modified sites"),
                       Line2D([], [], marker="o", mfc="white", mec=S.INK, lw=0, ms=4, label="background sites")],
              loc="upper left", bbox_to_anchor=(0.0, -0.2), ncol=2, borderaxespad=0)
    return fig, log


# ------------------------------------------------------------------------------------------------ Fig. 4 (old F2)
def rows_panel(ax, records, log, panel):
    y, ticks, labels = 0.0, [], []
    for verdict in F2.VERDICT_ORDER:
        block = sorted([r for r in records if r["verdict"] == verdict], key=F2.sort_key)
        if not block:
            continue
        # Verdict words live in the legend; leaving the row-group whitespace
        # unlabelled keeps the zero reference line and the claim labels clear.
        y -= 1.0
        for r in block:
            c = VC[verdict]
            if r["baseline_ci_low"] is not None:
                S.ci(ax, y + 0.24, r["baseline_log2_or"], r["baseline_ci_low"], r["baseline_ci_high"], S.INK, "o", False, 3.4)
            if r["matched_ci_low"] is not None:
                S.ci(ax, y, r["matched_log2_or"], r["matched_ci_low"], r["matched_ci_high"], c, "o", True, 3.8)
            if r["random_control_ci_low"] is not None:
                S.ci(ax, y - 0.24, r["random_control_log2_or"], r["random_control_ci_low"], r["random_control_ci_high"],
                     S.GREY, "x", True, 3.4)
            for k in ("baseline_log2_or", "baseline_ci_low", "baseline_ci_high", "matched_log2_or", "matched_ci_low",
                      "matched_ci_high", "random_control_log2_or", "random_control_ci_low", "random_control_ci_high"):
                if r[k] is not None:
                    log.add(panel, r["claim_id"], k, r[k + "__raw"], r[k], r["source_table"])
            ticks.append(y); labels.append(F2.label_of(r)); y -= 1.0
        y -= 0.3
    ax.axvline(0.0, color=S.INK, lw=S.LW, zorder=0)
    ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 0.5, 0.7)


def fig4():
    log = ValueLedger("Fig4")
    rows, no_baseline, _ = F2.collect()
    own = [r for r in rows if r["claim_id"] in F2.OWN_ESTIMAND]
    rest = [r for r in rows if r["claim_id"] not in F2.OWN_ESTIMAND]
    site = [r for r in rest if r["unit"] == "site"]
    prot = [r for r in rest if r["unit"] != "site"]
    fig = plt.figure(figsize=(width_mm * S.MM, 150 * S.MM))
    gs = fig.add_gridspec(3, 2, left=0.13, right=0.985, bottom=0.11, top=0.96, wspace=0.55, hspace=0.55,
                          width_ratios=[1, 1], height_ratios=[1, 0.26, 0.32])
    ax = fig.add_subplot(gs[:, 0]); rows_panel(ax, site, log, "a")
    ax.set_xlabel("log2 odds ratio"); ax.set_title("Site-level claims (%d)" % len(site)); letter(ax, "a", -0.02)
    ax = fig.add_subplot(gs[0, 1]); rows_panel(ax, prot, log, "b")
    ax.set_xlabel("log2 odds ratio"); ax.set_title("Protein-level claims (%d)" % len(prot)); letter(ax, "b", -0.02)
    ax = fig.add_subplot(gs[1, 1]); rows_panel(ax, own, log, "c")
    ax.set_xlabel("log2(observed / expected shared)"); ax.set_title("PERS-010, own estimand"); letter(ax, "c", -0.02)
    handles = [Line2D([], [], marker="o", color=S.INK, mfc="white", lw=S.LW, ms=3.6, label="baseline (paper's caliber)"),
               Line2D([], [], marker="o", color=S.INK, lw=S.LW, ms=3.8, label="after detectability matching"),
               Line2D([], [], marker="x", color=S.GREY, lw=S.LW, ms=3.6, label="same-size random control")]
    handles += [Line2D([], [], marker="o", color=VC[v], lw=0, ms=4, label=VT[v]) for v in F2.VERDICT_ORDER
                if any(r["verdict"] == v for r in rows)]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.56, 0.02), ncol=2, columnspacing=1.0,
               handletextpad=0.5)
    return fig, log, no_baseline


# ------------------------------------------------------------------------------------------------ Fig. 5 (old F3)
def fig5():
    log = ValueLedger("Fig5")
    aud = read(RES / "self_audit_public_cohorts_2026-09-20.csv")
    abl = read(RES / "public_refit_ablation_2026-09-21.csv")
    five = read(RES / "self_audit_public_five_proteases_2026-09-20.csv")
    cohorts = [("human_PXD044043", "Human", HUMAN), ("rice_PXD072089", "Rice", RICE)]
    fig = plt.figure(figsize=(width_mm * S.MM, 128 * S.MM))
    gs = fig.add_gridspec(2, 2, left=0.205, right=0.985, bottom=0.08, top=0.95, hspace=0.5, wspace=0.8)

    ax = fig.add_subplot(gs[0, 0]); y, t, l = 0, [], []
    for coh, cn, col in cohorts:
        for st, sn in (("label_positive", "positives"), ("label_negative", "negatives")):
            for key, mn, c in (("v2_chem_deployed", "chemistry", col), ("v2_full_benchmark", "full", S.GREY)):
                r = next(x for x in aud if x["cohort"] == coh and x["model"] == key and x["band"] == "distal_6_12"
                         and x["stratum"] == st)
                S.ci(ax, y, float(r["auc"]), float(r["ci_low"]), float(r["ci_high"]), c)
                for f in ("auc", "ci_low", "ci_high"):
                    log.add("a", f"{coh}|{key}|{st}", f, r[f], float(r[f]), "self_audit_public_cohorts_2026-09-20.csv")
                t.append(y); l.append(f"{cn} {sn}, {mn}"); y -= 1
        y -= 0.6
    ax.axvline(0.5, color=S.GREY, lw=S.LW, ls=(0, (2, 2)))
    ax.set_yticks(t); ax.set_yticklabels(l); ax.tick_params(axis="y", length=0); ax.set_ylim(y + 1.3, 0.7)
    ax.set_xlabel("AUC, score vs distal cleavage band"); ax.set_title("Within each label stratum"); letter(ax, "a", -0.02)

    ax = fig.add_subplot(gs[0, 1]); y, t, l = 0, [], []
    arms = [("control", "deployed columns"), ("drop_w710", "minus 2 windows"), ("drop_all_basic", "minus all basic"),
            ("label_permuted", "labels permuted")]
    for coh, cn, col in cohorts:
        for arm, an in arms:
            r = next(x for x in abl if x["cohort"] == coh and x["arm"] == arm)
            S.ci(ax, y, float(r["log2_or"]), float(r["ci_low"]), float(r["ci_high"]), S.GREY if arm == "label_permuted" else col)
            for f in ("log2_or", "ci_low", "ci_high"):
                log.add("b", f"{coh}|{arm}", f, r[f], float(r[f]), "public_refit_ablation_2026-09-21.csv")
            t.append(y); l.append(f"{cn}, {an}"); y -= 1
        y -= 0.6
    ax.axvline(0, color=S.INK, lw=S.LW)
    ax.set_yticks(t); ax.set_yticklabels(l); ax.tick_params(axis="y", length=0); ax.set_ylim(y + 1.3, 0.7)
    ax.set_xlabel("Top-100 log2 odds ratio"); ax.set_title("Deleting named columns"); letter(ax, "b", -0.02)

    ax = fig.add_subplot(gs[1, 0])
    rows5 = [("control", "deployed columns"), ("drop_w710", "minus 2 windows"), ("drop_all_basic", "minus all basic"),
             ("detect_only", "25 detectability features"), ("full", "all 1,107 columns")]
    for off, (coh, cn, col) in zip((-0.15, 0.15), cohorts):
        xs = []
        for i, (arm, an) in enumerate(rows5):
            r = next(x for x in abl if x["cohort"] == coh and x["arm"] == arm)
            xs.append(float(r["within_protein_auc"]))
            log.add("c", f"{coh}|{arm}", "within_protein_auc", r["within_protein_auc"], xs[-1], "public_refit_ablation_2026-09-21.csv")
        ax.plot(xs, [i + off for i in range(len(rows5))], "o", color=col, ms=4, label=cn)
    ax.set_yticks(range(len(rows5))); ax.set_yticklabels([a[1] for a in rows5]); ax.invert_yaxis()
    ax.tick_params(axis="y", length=0); ax.set_xlim(0.75, 0.92)
    ax.set_xlabel("Within-protein AUC"); ax.set_title("Detectability-only features"); letter(ax, "c", -0.02)
    ax.legend(loc="lower left", ncol=2)

    ax = fig.add_subplot(gs[1, 1])
    prot = ["Trypsin", "LysC", "GluC", "AspN", "Chymotrypsin"]
    for off, (coh, cn, col) in zip((-0.15, 0.15), cohorts):
        for i, p in enumerate(prot):
            r = next(x for x in five if x["cohort"] == coh and x["protease"] == p and x["model"] == "v2_chem_deployed"
                     and x["band"] == "distal_6_12" and x["stratum"] == "label_positive")
            S.ci(ax, i + off, float(r["auc"]), float(r["ci_low"]), float(r["ci_high"]), col)
            for f in ("auc", "ci_low", "ci_high"):
                log.add("d", f"{coh}|{p}", f, r[f], float(r[f]), "self_audit_public_five_proteases_2026-09-20.csv")
    ax.axvline(0.5, color=S.GREY, lw=S.LW, ls=(0, (2, 2)))
    ax.set_yticks(range(5)); ax.set_yticklabels(prot); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("AUC, label-positive stratum"); ax.set_title("Five cleavage rules"); letter(ax, "d", -0.02)
    return fig, log


def require_matplotlib_panel_alignment(fig):
    """Fail closed if any panel has invalid normalized geometry before export."""
    bad = []
    for i, ax in enumerate(fig.get_axes()):
        p = ax.get_position()
        if p.width <= 0 or p.height <= 0 or p.x0 < -1e-6 or p.y0 < -1e-6 or p.x1 > 1.000001 or p.y1 > 1.000001:
            bad.append((i, round(p.x0, 5), round(p.y0, 5), round(p.x1, 5), round(p.y1, 5)))
    if bad:
        raise RuntimeError(f"panel alignment failed: {bad}")


def gates_for(fig, name):
    # Render-time panel-alignment gate; keep this check in the plotting source.
    require_matplotlib_panel_alignment(fig)
    res = CHK.audit_figure(name, fig)
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); fw, fh = fig.canvas.get_width_height()
    drawn = list(fig.texts)
    for leg in fig.legends:
        drawn += list(leg.get_texts())
    for a in fig.get_axes():
        drawn += [a.title, a._left_title, a._right_title, a.xaxis.label, a.yaxis.label] + list(a.texts)
        if a.get_legend():
            drawn += list(a.get_legend().get_texts())
        for axis, lim in ((a.xaxis, a.get_xlim()), (a.yaxis, a.get_ylim())):
            lo, hi = min(lim), max(lim)
            for tk in axis.get_major_ticks():
                if lo - 1e-9 <= tk.get_loc() <= hi + 1e-9 and tk.label1.get_visible():
                    drawn.append(tk.label1)
    drawn = [t for t in drawn if t.get_visible() and (t.get_text() or "").strip()]
    outside = []
    for t in drawn:
        e = t.get_window_extent(rend)
        if e.x0 < -0.5 or e.y0 < -0.5 or e.x1 > fw + 0.5 or e.y1 > fh + 0.5:
            outside.append(t.get_text()[:40])
    short = [round(a.get_position().height * fig.get_figheight() * 25.4, 1) for a in fig.get_axes()
             if a.get_position().height * fig.get_figheight() * 25.4 < 20]
    return dict(collisions=res["n_collisions"], pairs=[(c["text_a"], c["text_b"]) for c in res["collisions"]],
                min_lw=S.min_linewidth(fig), font_sizes=S.font_sizes(fig), outside=sorted(set(outside)),
                axes_heights_mm=[round(a.get_position().height * fig.get_figheight() * 25.4, 1) for a in fig.get_axes()],
                size_mm=[round(fig.get_figwidth() * 25.4, 1), round(fig.get_figheight() * 25.4, 1)])


def main():
    OUT.mkdir(exist_ok=True); SD.mkdir(exist_ok=True)
    ctl = plt.figure(); ctl.text(.5, .5, "ctl a", fontsize=10); ctl.text(.5, .5, "ctl b", fontsize=10)
    pc = CHK.audit_figure("ctl", ctl)["n_collisions"]; plt.close(ctl)
    report, ok = {}, pc > 0
    extra = {}
    for name, fn in (("Fig2_search_space", fig2), ("Fig3_three_axes", fig3), ("Fig4_claim_retests", fig4),
                     ("Fig5_self_audit_public", fig5)):
        out = fn()
        fig, log = out[0], out[1]
        if len(out) > 2:
            extra["Fig4_no_baseline_claims"] = [(r["claim_id"], r["verdict"]) for r in out[2]]
        g = gates_for(fig, name)
        g["values_logged"] = len(log.rows)
        fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.svg")
        fig.savefig(OUT / f"{name}.png", dpi=600); fig.savefig(OUT / f"{name}.tiff", dpi=600)
        plt.close(fig)
        with open(SD / f"Source_Data_{name}.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(log.rows[0])); w.writeheader(); w.writerows(log.rows)
        g_ok = (g["collisions"] == 0 and g["min_lw"] >= 1.0 and set(g["font_sizes"]) <= {7.0, 8.0} and not g["outside"])
        g["pass"] = g_ok; ok = ok and g_ok
        report[name] = g
        print(name, json.dumps({k: v for k, v in g.items() if k != "pairs"}), g["pairs"][:4])
    (RES / "figs_2to5_nc_style_2026-09-21_audit.json").write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=hashlib.sha256(SELF).hexdigest(),
        imported=dict(f1=F1_SHA, f2=F2_SHA, style=S_SHA, detector=CHK_SHA), positive_control=pc,
        renumbering={"Fig2": "old F4", "Fig3": "old F1", "Fig4": "old F2", "Fig5": "old F3 (public)"},
        figures=report, extra=extra, all_pass=ok), indent=2, ensure_ascii=False))
    print("all_pass", ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
