"""Figure 3, revision 2026-09-30 (final): what a search can decide between sulfide and cysteine dioxidation (PXD015307).

Redrawn after revision item H_artifact5_docs (report section 7c, 'Figure 3'; Supplemental Note 6). Panels a, b and d
are POST HOC; panel c keeps its content and its design.

FINAL (results/figures_final/Fig3): results/figures_rev/Fig3/fig3_rev.py with four layout fixes found in verification;
no value, source or Source Data row changes (the Source Data is byte-identical to the candidate's):
  1  panel titles 7 pt regular (letters 8 pt bold), as in every other figure of the paper (Figs 1-7 as submitted);
     the candidate had 8-pt titles
  2  panel a: the '(%)' key sat on a white knock-out box across the dotted 1,776- and 3,946-Da lines; it now sits in
     the empty band between the 4.5- and 10-ppm lines, right-aligned with the other right-hand labels, with no box
  3  panel a strip: '7,523' ran over (knocked out) the dotted 8,879-Da line; the strip values are now 5.9 pt (as the
     tied/n counts of panel d), set 1 pt after the bar end, with no box
  4  panel d: the grey 'all' divider crossed both series labels ('...insulin files)' and '...Fig-1D)'); the axes are
     widened to end where panel b's axes end (6.36 in) and the 0.02-Da label starts at x = 0.8 instead of 1.2
A PDF-level gate (text ink boxes against every stroked segment, marker and bar) now fails the build on any overlap.

  a  the 0.017758-Da separation in ppm of peptide mass against +/-10, +/-4.5 and +/-2 ppm precursor windows
     (window criterion, unbiased measurement); below, each re-searched file's 2-SD mass range (accuracy criterion)
  b  the deposit's own insulin-file precursor errors (53 PSMs) coloured by insulin chain; below, the four
     sulfide-assigned PSMs against the errors expected under sulfide and under dioxidation, relative to B-chain
     and A-chain baselines
  c  PSMs at 1% FDR per file in the re-search; the +32-class count is zero in all five files (current design kept)
  d  share of single-cysteine spectra, with both candidates inside the +/-10-ppm window, whose two interpretations
     tie at the best xcorr, by fragment binning (0.02-Da FTMS MS2 vs 1.0005-Da ion-trap MS2) and best E-value

DATA CONTRACT. No statistic is recomputed. Step 1 copies every value that is drawn or printed, as the verbatim
string of a stored table, into a ledger that becomes Source_Data_Fig3_search_space.csv; step 2 draws only from
that ledger. One element is evaluated rather than copied: the separation curve of panel a, the exact identity
sep(M) = K / M ppm with K = t x M*(t) read from the stored window rows (all four stored rows give the same K,
asserted); its grid points are written to the Source Data too. Checks that compare stored tables with each other
(counts, readings, the unchanged panel c) are gates only and feed nothing into the drawing.

Baselines shown in panel b (from deposit_sulfide_psms_baseline_sensitivity.csv):
  A-chain sulfide PSMs: D1 'file, all non-sulfide PSMs (median, SD)' = the 12 B-chain PSMs of the CTH-with-Cys
                        file (the B-chain baseline), and D6b 'four files, same chain (file SD as stand-in)' = the
                        deposit's two A-chain PSMs (the A-chain baseline; the SD of two PSMs is not used, see README)
  B-chain sulfide PSMs: D6a 'four files, same chain (SD of those PSMs)' = the 47 B-chain PSMs (the B-chain baseline;
                        these two PSMs read the same under all eight stored baselines)
Filled expectation markers mean |z| < 2 (stored z_if_sulfide / z_if_dioxidation), the reading rule of the item.

Run:  python fig3_rev.py            -> writes into this script's folder
      python fig3_rev.py --out DIR  -> drafts elsewhere
Interpreter: system Python 3.10, matplotlib 3.10, PyMuPDF (fitz) for the PDF checks.
"""
import argparse
import csv
import hashlib
import json
import pathlib
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")

HERE = pathlib.Path(__file__).resolve().parent
WORK = pathlib.Path("C:/Users/admin/Desktop/小论文/_cys_repo_work")
H5 = WORK / "public/revision_2026-09-30/results/H_artifact5_docs"
REPO = WORK / "repo/results"
MCP_SD = pathlib.Path("C:/Users/admin/Desktop/小论文/巯基化/MCP/source_data")

SOURCES = {   # key -> (absolute path, label written into source_table)
    "decid": (H5 / "precursor_accuracy_decidability.csv", "H_artifact5_docs/precursor_accuracy_decidability.csv"),
    "sens": (H5 / "deposit_sulfide_psms_baseline_sensitivity.csv",
             "H_artifact5_docs/deposit_sulfide_psms_baseline_sensitivity.csv"),
    "nonsulf": (H5 / "deposit_nonsulfide_by_arm_chain_charge.csv",
                "H_artifact5_docs/deposit_nonsulfide_by_arm_chain_charge.csv"),
    "reass": (H5 / "deposit_sulfide_psms_reassessed.csv", "H_artifact5_docs/deposit_sulfide_psms_reassessed.csv"),
    "robust": (H5 / "deposit_sulfide_psms_robust_reading.csv",
               "H_artifact5_docs/deposit_sulfide_psms_robust_reading.csv"),
    "summ32": (H5 / "deposit_plus32_summary.json", "H_artifact5_docs/deposit_plus32_summary.json"),
    "fragarith": (H5 / "fragment_bin_arithmetic.csv", "H_artifact5_docs/fragment_bin_arithmetic.csv"),
    "ties": (H5 / "single_cys_tie_rate_when_both_scored.csv",
             "H_artifact5_docs/single_cys_tie_rate_when_both_scored.csv"),
    "psms": (REPO / "search_space_mass_accuracy.csv", "repo/results/search_space_mass_accuracy.csv"),
    "summary": (REPO / "pxd015307_research_summary.csv", "repo/results/pxd015307_research_summary.csv"),
    "curve": (None, "derived: K/M ppm, K = tolerance x mass_limit_da of the stored window rows "
                    "(H_artifact5_docs/precursor_accuracy_decidability.csv)"),
}
CURRENT_SD = MCP_SD / "Source_Data_Fig3_search_space.csv"      # read-only: panel c must be unchanged

# ------------------------------------------------------------------------------------------------ house style
W_IN = 7.2                     # exactly 518.4 pt
FS, FS_S, FS_T = 7.0, 6.5, 8.0
LW_AX, LW = 0.8, 1.0
INK, AXIS, MUTE, GREY, REF, BAND = "#1f2326", "#4a4f55", "#6f757c", "#9aa1a8", "#c3c8ce", "#f1f3f5"
BLUE, ORANGE, VERMIL = "#1f6aa5", "#dc9a2b", "#cc5a36"      # v2 palette (Okabe-Ito hues, deepened)
CHAIN = {"B": BLUE, "A": ORANGE}
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"], "font.size": FS,
    "axes.titlesize": FS_T, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
    "legend.fontsize": FS_S, "axes.linewidth": LW_AX, "lines.linewidth": LW,
    "xtick.major.width": LW_AX, "ytick.major.width": LW_AX, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0, "axes.labelpad": 3.0,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": AXIS,
    "xtick.color": AXIS, "ytick.color": AXIS, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "text.color": INK, "axes.labelcolor": INK, "axes.unicode_minus": True, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular",
})
ARM_LABEL = {"Fig-1D": "Fig-1D", "Fig4-L_CTH_Cys": "CTH Cys", "Fig4-L_NaHS": "NaHS", "Fig4-L_control": "control",
             "Fig4-L_heat_inactive_CTH": "heat inactive CTH"}     # the labels of the current panel c
DEP_ARMS = ["CTH_Cys", "NaHS", "control", "heat_inactive_CTH"]     # deposit arm names, same order as panel c


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def read_csv(key):
    with open(SOURCES[key][0], encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fmt_int(x):
    return f"{x:,.0f}"


# ------------------------------------------------------------------------------------------------ step 1: ledger
class Ledger:
    """Every drawn or printed value, as the verbatim string of its stored table."""

    def __init__(self):
        self.rows, self.index, self.used = [], {}, set()

    def add(self, panel, row, field, raw, key, spec=""):
        k = (panel, row, field)
        assert k not in self.index, k
        self.index[k] = len(self.rows)
        self.rows.append(dict(figure="Fig3", panel=panel, row=row, field=field, value=str(raw),
                              source_table=SOURCES[key][1], specification_label=spec))

    def v(self, panel, row, field):
        k = (panel, row, field)
        self.used.add(k)
        return float(self.rows[self.index[k]]["value"])

    def s(self, panel, row, field):
        k = (panel, row, field)
        self.used.add(k)
        return self.rows[self.index[k]]["value"]


WIN_CRIT = "window: alternative outside the precursor window (M < M*)"
ACC_CRIT = "accuracy: separation >= 4 SD (2 SD each side of the midpoint)"
SHOW_B = {  # (chain of the sulfide PSM, baseline_definition) -> short label printed in panel b
    ("A", "D1 file, all non-sulfide PSMs (median, SD)"): "B chain, same file",
    ("A", "D6b four files, same chain (file SD as stand-in)"): "A chain",
    ("B", "D6a four files, same chain (SD of those PSMs)"): "B chain",
}
BASE_COLOUR = {"B chain, same file": BLUE, "A chain": ORANGE, "B chain": BLUE}


def collect():
    L = Ledger()
    gates = {}

    # ---- panel a: window criterion (unbiased) and the per-file 2-SD mass range
    dec = read_csv("decid")
    win = [r for r in dec if r["criterion"] == WIN_CRIT]
    ks = [float(r["parameter"].split()[1]) * float(r["mass_limit_da"]) for r in win]
    assert len(win) == 4 and (max(ks) - min(ks)) / min(ks) < 1e-12, ks
    gates["curve_constant_K_ppm_Da_all_window_rows_agree"] = [min(ks), max(ks)]
    for t in ("10.0", "4.5", "2.0"):
        r = [x for x in win if x["parameter"] == f"tolerance {t} ppm"]
        assert len(r) == 1
        r = r[0]
        row = f"window +/-{t} ppm"
        spec = f"{r['criterion']}; {r['context']}"
        L.add("a", row, "tolerance_ppm", t, "decid", spec + " (value from the 'parameter' cell)")
        L.add("a", row, "mass_limit_da", r["mass_limit_da"], "decid", spec)
        L.add("a", row, "fraction_theoretical_cys_peptides_at_or_below_plus32_only",
              r["fraction_theoretical_cys_peptides_at_or_below_plus32_only"], "decid", spec)
    r10 = [x for x in win if x["parameter"] == "tolerance 10.0 ppm"][0]
    k = 10.0 * float(r10["mass_limit_da"])
    fa = [r for r in read_csv("fragarith") if r["fragment_charge"] == "1"]
    assert len({r["shift_mz"] for r in fa}) == 1 and abs(float(fa[0]["shift_mz"]) - k * 1e-6) < 1e-7
    L.add("a", "separation label", "separation_da", fa[0]["shift_mz"], "fragarith",
          "printed curve label; dioxidation minus sulfide (shift_mz at fragment charge 1)")
    grid, m = [], 1100
    while m <= 13000:                              # 1.1-13 kDa; finer where the curve bends
        grid.append(m)
        m += 20 if m < 2000 else 50 if m < 5000 else 100
    for i, m in enumerate(grid, 1):
        row = f"separation curve point {i:03d}"
        L.add("a", row, "mass_da", str(m), "curve", "grid point of the drawn curve")
        L.add("a", row, "separation_ppm", repr(k / m), "curve",
              f"K / mass_da, K = 10.0 x {r10['mass_limit_da']} = {k!r} ppm x Da")
    acc = [r for r in dec if r["criterion"] == ACC_CRIT and r["context"].startswith("PXD015307 re-search ")]
    summ = read_csv("summary")
    arms = [r["arm"] for r in summ]
    assert len(acc) == 5
    for arm in arms:
        r = [x for x in acc if x["context"].startswith(f"PXD015307 re-search {arm} (")]
        assert len(r) == 1, arm
        r = r[0]
        L.add("a", f"2-SD mass range|{arm}", "mass_limit_da", r["mass_limit_da"], "decid",
              f"{r['criterion']}; {r['context']}; {r['parameter']}")

    # ---- panel b: the deposit's 53 insulin-file PSMs, and the four sulfide PSMs against two baselines
    psms = read_csv("psms")
    reass = read_csv("reass")
    seq2chain = {r["sequence"]: r["chain"] for r in reass}
    assert set(seq2chain.values()) == {"A", "B"} and len(seq2chain) == 2
    assert len(psms) == 53 and all(p["sequence"] in seq2chain for p in psms)
    js = json.loads(SOURCES["summ32"][0].read_text(encoding="utf-8"))
    L.add("b", "deposit PSMs|control", "n_psms", js["targeted_arms"]["control"]["n_psms"], "summ32",
          "targeted_arms.control.n_psms (the control consensus file carries no PSM)")
    for i, p in enumerate(psms, 1):
        ch = seq2chain[p["sequence"]]
        row = f"deposit PSM {i:02d}|{p['arm']}|{ch} chain|z={p['charge']}|{p['placed']}"
        L.add("b", row, "observed_ppm", p["observed_ppm"], "psms",
              "the deposit's TargetPsms.DeltaMassInPPM (stored extraction); chain from the sequence "
              "(H_artifact5_docs/deposit_sulfide_psms_reassessed.csv); "
              + ("sulfide-assigned" if p["has_sulfide"] == "True" else "no sulfide"))
    # gate: non-sulfide counts per arm x chain x charge equal the stored group table
    grp = {}
    for p in psms:
        if p["has_sulfide"] == "False":
            kk = (p["arm"], seq2chain[p["sequence"]], p["charge"])
            grp[kk] = grp.get(kk, 0) + 1
    ns = {(r["arm"], r["chain"], r["charge"]): int(r["n"]) for r in read_csv("nonsulf")}
    assert grp == ns, (grp, ns)
    gates["nonsulfide_counts_match_deposit_nonsulfide_by_arm_chain_charge"] = True
    sens = read_csv("sens")
    shown = [r for r in sens if (r["chain"], r["baseline_definition"]) in SHOW_B]
    assert len(shown) == 6, len(shown)
    sulf = [p for p in psms if p["has_sulfide"] == "True"]
    for r in shown:
        assert any(p["arm"] == r["arm"] and p["observed_ppm"] == r["observed_ppm"] for p in sulf), r
        zs, zo = float(r["z_if_sulfide"]), float(r["z_if_dioxidation"])
        rule = ("sulfide only" if abs(zs) < 2 <= abs(zo) else "either" if abs(zs) < 2 and abs(zo) < 2
                else "dioxidation only" if abs(zo) < 2 else "neither")
        assert rule == r["reading"], (r, rule)          # the stored reading follows the stored z (gate only)
        row = (f"sulfide PSM|{r['arm']}|{r['chain']} chain|z={r['charge']}|{r['placed']}|"
               f"{r['observed_ppm']}|baseline {r['baseline_definition'].split()[0]}")
        spec = f"{r['baseline_definition']}; baseline PSMs {r['baseline_chains']}; spread from: {r['spread_source']}"
        for f in ("observed_ppm", "n_baseline_psms", "centre_ppm", "expected_ppm_if_dioxidised_but_reported_as_sulfide",
                  "z_if_sulfide", "z_if_dioxidation", "reading"):
            L.add("b", row, f, r[f], "sens", spec)
    gates["panel_b_readings_follow_stored_z_rule"] = True
    rob = {f"{r['arm']} {r['chain']} chain {r['observed_ppm']} ppm": r["robust_summary"] for r in read_csv("robust")}
    assert len(rob) == 4
    gates["panel_b_robust_summaries"] = rob

    # ---- panel c: unchanged content
    for r in summ:
        L.add("c", r["arm"], "n_psms_at_1pct_fdr", r["n_psms_at_1pct_fdr"], "summary")
        L.add("c", r["arm"], "n_plus32_cys_psms", r["n_plus32_cys_psms"], "summary")
    with open(CURRENT_SD, encoding="utf-8-sig", newline="") as fh:
        cur = {(r["row"], r["field"]): r["value"] for r in csv.DictReader(fh) if r["panel"] == "c"}
    new = {(x["row"], x["field"]): x["value"] for x in L.rows if x["panel"] == "c"}
    assert cur == new, (cur, new)
    gates["panel_c_values_identical_to_current_source_data"] = True

    # ---- panel d: tie share among single-Cys spectra with both candidates in the window
    ties = read_csv("ties")
    assert len(ties) == 10
    for r in ties:
        assert abs(float(r["share_tied"]) - int(r["n_tied"]) / int(r["n_both_scored"])) < 1e-12
        row = f"{r['group']}|best E-value {r['evalue_bin']}"
        for f in ("n_both_scored", "n_tied", "share_tied"):
            L.add("d", row, f, r[f], "ties", "single-cysteine spectra whose best interpretation carries one +32-class "
                  "modification and whose alternative lies inside the +/-10-ppm window (both scored)")
    return L, gates


# ------------------------------------------------------------------------------------------------ helpers
def swarm(xs_pt, diam_pt, max_off_pt):
    """Deterministic beeswarm (layout only, no data change): points given in points, placed in x order (larger
    markers first at equal x); each takes the smallest |vertical offset| at which it overlaps no placed marker.
    Returns the offsets and the worst overlap as a fraction of the touching distance (0 = no overlap)."""
    placed = []
    order = sorted(range(len(xs_pt)), key=lambda i: (-diam_pt[i], xs_pt[i], i))
    res = [0.0] * len(xs_pt)
    step = 0.4
    cands = [0.0]
    j = 1
    while j * step <= max_off_pt + 1e-9:
        cands += [j * step, -j * step]
        j += 1

    def clash(x, c, d):
        return max([0.0] + [1 - ((x - px) ** 2 + (c - py) ** 2) ** 0.5 / (0.5 * (d + pd)) for px, py, pd in placed])

    for i in order:
        x, d = xs_pt[i], diam_pt[i]
        best = next((c for c in cands if clash(x, c, d) <= 0.0), None)
        if best is None:        # no free slot inside the band: take the least-overlapping one
            best = min(cands, key=lambda c: (clash(x, c, d), abs(c)))
        placed.append((x, best, d))
        res[i] = best
    worst = 0.0
    for a in range(len(placed)):
        for b in range(a + 1, len(placed)):
            (xa, ya, da), (xb, yb, db) = placed[a], placed[b]
            worst = max(worst, 1 - ((xa - xb) ** 2 + (ya - yb) ** 2) ** 0.5 / (0.5 * (da + db)))
    return res, worst


def clean_y(ax):
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["left"].set_visible(False)


def label_panels(fig, panels, gap_mm=1.4, above_mm=2.4):
    """(axes group, letter, title, top axes): letter at the group's leftmost tight-box edge (a column shares one x),
    title after it on the same baseline, both above the top axes."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    fw, fh = fig.get_figwidth() * 25.4, fig.get_figheight() * 25.4
    lefts = []
    for axes, letter, title, top in panels:
        x0 = min(a.get_tightbbox(r).transformed(inv).x0 for a in axes)
        lefts.append((x0, min(a.get_position().x0 for a in axes)))
    out = []
    for (axes, letter, title, top), (x0, fx0) in zip(panels, lefts):
        col = [lx for lx, lf in lefts if abs(lf - fx0) * fw < 25]
        x = max(min(col), 0.004)
        y = top.get_position().y1 + above_mm / fh
        t = fig.text(x, y, letter, fontsize=FS_T, fontweight="bold", ha="left", va="baseline", color=INK)
        fig.canvas.draw()
        w = t.get_window_extent(r).transformed(inv).width
        out.append(fig.text(x + w + gap_mm / fw, y, title, fontsize=FS, ha="left", va="baseline", color=INK))
    return out


# ------------------------------------------------------------------------------------------------ step 2: draw
def draw(L):
    H_IN = 5.02
    fig = plt.figure(figsize=(W_IN, H_IN))

    def ax_in(x0, y0, w, h):      # inches from the bottom-left corner
        return fig.add_axes([x0 / W_IN, y0 / H_IN, w / W_IN, h / H_IN])

    # geometry (inches)
    AX_L, AX_W = 0.98, 2.30                  # panel a (curve and strip share x); panel c starts at AX_L too
    B_L, B_W = 4.56, 1.80                    # panel b
    ROW1_TOP = H_IN - 0.30
    A_CURVE_H, A_GAP, A_STRIP_H = 1.46, 0.20, 0.74
    ROW1_BOT = ROW1_TOP - (A_CURVE_H + A_GAP + A_STRIP_H)
    ROW2_H = 1.08
    ROW2_TOP = ROW1_BOT - 0.72

    # ================================================================ a
    ax = ax_in(AX_L, ROW1_TOP - A_CURVE_H, AX_W, A_CURVE_H)
    axs = ax_in(AX_L, ROW1_BOT, AX_W, A_STRIP_H)
    XMAX, YMAX = 13000, 16.0
    pts = sorted((L.v("a", r["row"], "mass_da"), L.v("a", r["row"], "separation_ppm"))
                 for r in L.rows if r["panel"] == "a" and r["field"] == "mass_da")
    mx, sp = [p[0] for p in pts], [p[1] for p in pts]
    ax.fill_between(mx, sp, YMAX, color=VERMIL, alpha=0.08, lw=0, zorder=0)
    ax.plot(mx, sp, color=VERMIL, lw=1.4, zorder=3, solid_capstyle="butt")
    ax.text(1560, 13.2, f"{L.s('a', 'separation label', 'separation_da')} Da in ppm of M", color=VERMIL,
            fontsize=FS_S, ha="left", va="center")
    ax.text(12850, 14.6, "both candidates inside the window", color=MUTE, fontsize=FS_S, ha="right",
            va="center")
    for t in ("10.0", "4.5", "2.0"):
        row = f"window +/-{t} ppm"
        tol, mstar = L.v("a", row, "tolerance_ppm"), L.v("a", row, "mass_limit_da")
        frac = L.v("a", row, "fraction_theoretical_cys_peptides_at_or_below_plus32_only")
        ax.axhline(tol, color=GREY, lw=LW_AX, ls=(0, (3, 2)), zorder=1)
        ax.plot([mstar, mstar], [0, tol], color=GREY, lw=LW_AX, ls=(0, (1, 1.6)), zorder=1)
        axs.axvline(mstar, color=GREY, lw=LW_AX, ls=(0, (1, 1.6)), zorder=1)
        ax.plot([mstar], [tol], "o", ms=3.4, color=INK, mew=0, zorder=5)
        pct = ">99.9%" if 0.999 < frac < 1.0 else f"{frac:.0%}"
        lab = f"±{t.rstrip('0').rstrip('.')} ppm: {fmt_int(mstar)} Da ({pct})"
        if t == "2.0":      # no room right of the last crossing: right-aligned, between the 2- and 4.5-ppm lines
            ax.text(12850, 3.35, lab, color=INK, fontsize=FS_S, ha="right", va="center", zorder=6)
        else:
            ax.text(mstar + 200, tol + 0.3, lab, color=INK, fontsize=FS_S, ha="left", va="bottom", zorder=6)
    # key to the percentages: in the empty band between the 4.5- and 10-ppm lines, right of the curve and of every
    # dotted M* line (the bottom-left corner it used to occupy is crossed by the 1,776- and 3,946-Da lines)
    ax.text(12850, 7.6, "(%): theoretical Cys peptides ≤ that mass", color=MUTE, fontsize=5.9,
            ha="right", va="center", zorder=4)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, YMAX)
    ax.set_yticks([0, 5, 10, 15])
    ax.set_ylabel("Separation (ppm)")
    ax.set_xticks(range(0, 12001, 2000))
    ax.tick_params(axis="x", labelbottom=False)
    # strip: mass up to which each file's precursor precision separates the candidates by 2 SD on each side
    arms = [r["row"] for r in L.rows if r["panel"] == "c" and r["field"] == "n_psms_at_1pct_fdr"]
    for i, arm in enumerate(arms):
        m2 = L.v("a", f"2-SD mass range|{arm}", "mass_limit_da")
        axs.barh(i, m2, height=0.60, color=GREY, lw=0, zorder=2)
        # 5.9 pt, 1 pt after the bar end, no knock-out box: '7,523' then ends clear of the dotted 8,879-Da line
        axs.annotate(fmt_int(m2), xy=(m2, i), xytext=(1.0, 0), textcoords="offset points", color=MUTE,
                     fontsize=5.9, va="center", ha="left", zorder=4)
    axs.set_yticks(range(len(arms)))
    axs.set_yticklabels([ARM_LABEL[a] for a in arms])
    axs.set_ylim(len(arms) - 0.45, -0.55)
    clean_y(axs)
    axs.text(0.0, 1.04, "mass range separated at 2 SD (each file’s own precision)", color=MUTE,
             fontsize=FS_S, ha="left", va="bottom", transform=axs.transAxes)
    axs.set_xlim(0, XMAX)
    axs.set_xticks(range(0, 12001, 2000))
    axs.set_xticklabels([fmt_int(x) for x in range(0, 12001, 2000)])
    axs.set_xlabel("Peptide neutral mass (Da)")

    # ================================================================ b
    axb = ax_in(B_L, ROW1_BOT, B_W, ROW1_TOP - ROW1_BOT)
    XB0, XB1 = -8.5, 14.5
    y_arm = {"CTH_Cys": 0.0, "NaHS": 1.2, "control": 2.2, "heat_inactive_CTH": 3.0}
    srows = [r["row"] for r in L.rows if r["panel"] == "b" and r["row"].startswith("sulfide PSM|")
             and r["field"] == "observed_ppm"]
    psm_keys = []
    for rk in srows:
        key = "|".join(rk.split("|")[:6])
        if key not in psm_keys:
            psm_keys.append(key)
    assert len(psm_keys) == 4
    lower, spans, y = [], [], 4.3
    for key in psm_keys:
        subs = sorted([rk for rk in srows if rk.startswith(key + "|")], key=lambda s: s.endswith("D6b"))
        y_first = y
        for rk in subs:                                   # B-chain baseline first, then the A-chain baseline
            lower.append((y, rk))
            y += 0.6
        spans.append((y_first - 0.3, y - 0.3))
        y += 0.3
    y_end = y - 0.3
    for g, (s0, s1) in enumerate(spans):                 # pale band behind every other PSM (house row grouping)
        if g % 2 == 0:
            axb.add_patch(Rectangle((0, s0), 1, s1 - s0, transform=axb.get_yaxis_transform(), facecolor=BAND,
                                    edgecolor="none", lw=0, zorder=0))
    YB0 = -0.55
    swarm_overlap = {}
    axb.set_xlim(XB0, XB1)
    axb.set_ylim(y_end - 0.1, YB0)
    fig.canvas.draw()
    bb = axb.get_window_extent()
    ppx = bb.width / (XB1 - XB0) * 72 / fig.dpi                  # points per ppm
    ppy = bb.height / (y_end - 0.1 - YB0) * 72 / fig.dpi          # points per row unit
    for arm in DEP_ARMS:
        rows = [r["row"] for r in L.rows if r["panel"] == "b" and r["row"].startswith("deposit PSM ")
                and r["row"].split("|")[1] == arm]
        yy = y_arm[arm]
        if not rows:
            assert L.v("b", "deposit PSMs|control", "n_psms") == 0
            axb.text(XB0 + 0.3, yy, "no PSM in the deposit", color=MUTE, fontsize=FS_S, va="center", ha="left")
            continue
        xs = [L.v("b", rk, "observed_ppm") for rk in rows]
        sulf = [rk.endswith("Sulfide") for rk in rows]
        offs, worst = swarm([x * ppx for x in xs], [4.6 if s else 3.0 for s in sulf], max_off_pt=0.56 * ppy)
        swarm_overlap[arm] = round(worst, 3)
        for x, o, rk, s in zip(xs, offs, rows, sulf):
            ch = rk.split("|")[2][0]
            if s:
                axb.plot([x], [yy + o / ppy], "D", ms=4.1, mfc=CHAIN[ch], mec=INK, mew=0.8, zorder=5)
            else:
                axb.plot([x], [yy + o / ppy], "o", ms=3.0, mfc=CHAIN[ch], mec="white", mew=0.3, zorder=4)
    axb.axhline(3.6, color=REF, lw=LW_AX, zorder=0)
    axb.text(XB0 - 0.35, 3.9, "baseline (n):", color=MUTE, fontsize=FS_S, ha="right", va="center")
    for y, rk in lower:
        chain = rk.split("|")[2][0]
        bdef = [r["specification_label"].split(";")[0] for r in L.rows
                if r["row"] == rk and r["field"] == "centre_ppm"][0]
        lab = SHOW_B[(chain, bdef)]
        col = BASE_COLOUR[lab]
        obs = L.v("b", rk, "observed_ppm")
        es, eo = L.v("b", rk, "centre_ppm"), L.v("b", rk, "expected_ppm_if_dioxidised_but_reported_as_sulfide")
        zs, zo = L.v("b", rk, "z_if_sulfide"), L.v("b", rk, "z_if_dioxidation")
        axb.plot([es, eo], [y, y], color=col, lw=LW, alpha=0.5, zorder=2, solid_capstyle="butt")
        axb.plot([es], [y], "o", ms=4.0, mfc=col if abs(zs) < 2 else "white", mec=col, mew=1.0, zorder=3)
        axb.plot([eo], [y], "s", ms=3.6, mfc=col if abs(zo) < 2 else "white", mec=col, mew=1.0, zorder=3)
        axb.plot([obs], [y], "D", ms=4.1, mfc=CHAIN[chain], mec=INK, mew=0.8, zorder=5)
        n = L.v("b", rk, "n_baseline_psms")
        axb.text(XB0 - 0.35, y, f"{lab} ({n:.0f})", color=col, fontsize=FS_S, ha="right", va="center")
        axb.text(XB1 + 0.25, y, L.s("b", rk, "reading"), color=MUTE, fontsize=FS_S, ha="left", va="center")
    axb.set_yticks([y_arm[a] for a in DEP_ARMS])
    axb.set_yticklabels([ARM_LABEL["Fig4-L_" + a] for a in DEP_ARMS])
    clean_y(axb)
    axb.set_xticks([-5, 0, 5, 10])
    axb.set_xlabel("Precursor error in the deposit’s search (ppm)")
    hb = [Line2D([], [], ls="", marker="o", ms=3.0, mfc=BLUE, mec="white", mew=0.3, label="B chain"),
          Line2D([], [], ls="", marker="o", ms=3.0, mfc=ORANGE, mec="white", mew=0.3, label="A chain"),
          Line2D([], [], ls="", marker="D", ms=4.1, mfc="white", mec=INK, mew=0.8, label="sulfide-assigned"),
          Line2D([], [], ls="", marker="o", ms=4.0, mfc="white", mec=MUTE, mew=1.0, label="expected if sulfide"),
          Line2D([], [], ls="", marker="s", ms=3.6, mfc="white", mec=MUTE, mew=1.0,
                 label="expected if dioxidation"),
          Line2D([], [], ls="", marker="s", ms=3.6, mfc=MUTE, mec=MUTE, mew=1.0, label="filled: PSM within 2 SD")]
    axb.legend(handles=hb, loc="upper left", bbox_to_anchor=((9.2 - XB0) / (XB1 - XB0), 1.0), ncol=1,
               fontsize=FS_S, handletextpad=0.35, borderaxespad=0.0, labelspacing=0.42, handlelength=1.0)

    # ================================================================ c (current design)
    axc = ax_in(AX_L, ROW2_TOP - ROW2_H, 2.05, ROW2_H)
    names = [ARM_LABEL[a] for a in arms]
    psms_n = [L.v("c", a, "n_psms_at_1pct_fdr") for a in arms]
    plus = [L.v("c", a, "n_plus32_cys_psms") for a in arms]
    axc.barh(range(len(arms)), psms_n, color=GREY, height=0.62, lw=0, zorder=2)
    for i, (n, p) in enumerate(zip(psms_n, plus)):
        axc.text(n + 6, i, f"{n:.0f}", va="center", color=MUTE)
        axc.text(1.02, i, f"{p:.0f}", va="center", ha="center", color=VERMIL, transform=axc.get_yaxis_transform())
    axc.text(1.02, -0.48, "+32", ha="center", va="bottom", color=VERMIL, transform=axc.get_yaxis_transform())
    axc.set_yticks(range(len(arms)))
    axc.set_yticklabels(names)
    axc.invert_yaxis()
    clean_y(axc)
    axc.set_xlim(0, 300)
    axc.set_xticks([0, 100, 200, 300])
    axc.set_xlabel("PSMs at 1% FDR")

    # ================================================================ d
    axd = ax_in(B_L - 0.28, ROW2_TOP - ROW2_H, B_W + 0.28, ROW2_H)     # right edge = panel b's (6.36 in)
    bins = ["(0.0, 1.0]", "(1.0, 10.0]", "(10.0, 100.0]", "(100.0, 1000.1]"]
    blab = ["≤1", "1–10", "10–100", ">100"]
    groups = [("insulin files (ion-trap MS2, 1.0005-Da bins)", VERMIL, "s", "1.0005-Da bins (ion trap; insulin files)"),
              ("Fig-1D (FTMS MS2, 0.02-Da bins)", INK, "o", "0.02-Da bins (FTMS; Fig-1D)")]
    XALL = 4.8
    for gi, (g, col, mk, lab) in enumerate(groups):
        ys = [L.v("d", f"{g}|best E-value {b}", "share_tied") for b in bins]
        axd.plot(range(4), ys, color=col, lw=LW, zorder=2)
        axd.plot(range(4), ys, mk, color=col, ms=3.6, mew=0, zorder=3)
        ya = L.v("d", f"{g}|best E-value all", "share_tied")
        axd.plot([XALL], [ya], mk, color=col, ms=5.0, mew=0, zorder=3)
        nt, nb = L.v("d", f"{g}|best E-value all", "n_tied"), L.v("d", f"{g}|best E-value all", "n_both_scored")
        axd.text(XALL + 0.3, ya, f"{ya:.0%}", color=col, fontsize=FS_S, va="center")
        for xi, b in enumerate(bins + ["all"]):
            rk = f"{g}|best E-value {b}"
            n_b, n_t = L.v("d", rk, "n_both_scored"), L.v("d", rk, "n_tied")
            axd.annotate(f"{n_t:.0f}/{n_b:.0f}", xy=(XALL if b == "all" else xi, 0),
                         xycoords=axd.get_xaxis_transform(), xytext=(0, -12.5 - 7.0 * gi),
                         textcoords="offset points", color=col, fontsize=5.9, ha="center", va="top")
    axd.annotate("tied/n", xy=(-0.45, 0), xycoords=axd.get_xaxis_transform(), xytext=(-3, -12.5),
                 textcoords="offset points", color=MUTE, fontsize=5.9, ha="right", va="top")
    axd.text(-0.3, 0.70, groups[0][3], color=VERMIL, fontsize=FS_S, ha="left", va="center")
    axd.text(0.8, 0.08, groups[1][3], color=INK, fontsize=FS_S, ha="left", va="center")
    axd.axvline(4.2, color=REF, lw=LW_AX, zorder=0)
    axd.set_xlim(-0.45, XALL + 0.45)
    axd.set_ylim(0, 1.05)
    axd.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axd.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    axd.set_xticks(list(range(4)) + [XALL])
    axd.set_xticklabels(blab + ["all"])
    axd.set_xlabel("Best E-value of the spectrum", labelpad=17)
    axd.set_ylabel("Share tied at best score")

    titles = label_panels(fig, [
        ([ax, axs], "a", "Separation against the precursor window", ax),
        ([axb], "b", "Deposit sulfide PSMs: the reading depends on the baseline", axb),
        ([axc], "c", "No +32 identification", axc),
        ([axd], "d", "Score ties follow the fragment binning", axd)])
    return fig, dict(a=ax, a_strip=axs, b=axb, c=axc, d=axd, swarm_overlap=swarm_overlap), titles


# ------------------------------------------------------------------------------------------------ gates
def text_gates(fig):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, Hh = fig.bbox.width, fig.bbox.height
    texts = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = sorted({round(t.get_fontsize(), 2) for t in texts})
    fams = sorted({t.get_fontname() for t in texts})
    boxes = [(t, t.get_window_extent(r)) for t in texts]
    outside = [t.get_text() for t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > Hh + 0.5]
    coll = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.x0 < b.x1 - 1 and b.x0 < a.x1 - 1 and a.y0 < b.y1 - 1 and b.y0 < a.y1 - 1:
                coll.append((boxes[i][0].get_text(), boxes[j][0].get_text()))
    return dict(font_sizes=sizes, font_families=fams, text_outside_canvas=outside, text_collisions=coll)


def pdf_gates(pdf):
    import fitz
    d = fitz.open(pdf)
    p = d[0]
    fonts = sorted({f[3] for f in d.get_page_fonts(0)})
    spans = [s for b in p.get_text("dict")["blocks"] for ln in b.get("lines", []) for s in ln["spans"]
             if s["text"].strip()]
    out = dict(page_pt=[round(p.rect.width, 3), round(p.rect.height, 3)], fonts=fonts,
               min_span_size_pt=round(min(s["size"] for s in spans), 3),
               span_fonts=sorted({s["font"] for s in spans}), n_pages=len(d))
    d.close()
    return out


def pdf_overlaps(pdf):
    """Text ink boxes (per word) of the written PDF against each other and against every drawn mark: stroked
    segments with half their line width, filled markers and bars with their box. Horizontal words use Arial metrics
    (cap height 0.716 em above the baseline, descender 0.212 em below it only for words with descenders); rotated
    words use their glyph boxes. Large fills (axes backgrounds, the shaded region of panel a, the row bands of
    panel b) are backgrounds and are not counted."""
    import math
    import fitz
    d = fitz.open(pdf)
    p = d[0]
    words = []
    for b in p.get_text("rawdict")["blocks"]:
        for ln in b.get("lines", []):
            horiz = abs(ln["dir"][0] - 1) < 1e-6
            for s in ln["spans"]:
                cur = []
                for c in s["chars"] + [None]:
                    if c is not None and c["c"].strip():
                        cur.append(c)
                        continue
                    if cur:
                        x0, x1 = min(c_["bbox"][0] for c_ in cur), max(c_["bbox"][2] for c_ in cur)
                        txt = "".join(c_["c"] for c_ in cur)
                        if horiz:
                            oy, sz = s["origin"][1], s["size"]
                            top = oy - (0.74 if any(ch in "bdfhklt()[]{}|/%" for ch in txt) else 0.716) * sz
                            bot = oy + (0.212 * sz if any(ch in "gjpqy(),;[]{}|/_Q" for ch in txt) else 0.0)
                            r = fitz.Rect(x0, top, x1, bot)
                        else:
                            r = fitz.Rect(x0, min(c_["bbox"][1] for c_ in cur), x1, max(c_["bbox"][3] for c_ in cur))
                        words.append((txt, r))
                    cur = []

    def seg_dist(a, b, r):
        if fitz.Rect(min(a.x, b.x), min(a.y, b.y), max(a.x, b.x) + 1e-9, max(a.y, b.y) + 1e-9).intersects(r):
            n = 64                                        # sampled test for a crossing segment
            for i in range(n + 1):
                x, y = a.x + (b.x - a.x) * i / n, a.y + (b.y - a.y) * i / n
                if r.x0 <= x <= r.x1 and r.y0 <= y <= r.y1:
                    return 0.0
        pts = [(a.x + (b.x - a.x) * i / 32, a.y + (b.y - a.y) * i / 32) for i in range(33)]
        return min(math.hypot(max(r.x0 - x, 0, x - r.x1), max(r.y0 - y, 0, y - r.y1)) for x, y in pts)

    tt, tm = [], []
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            inter = fitz.Rect(words[i][1]) & words[j][1]
            if not inter.is_empty and inter.width > 0.2 and inter.height > 0.2:
                tt.append((words[i][0], words[j][0]))
    for dr in p.get_drawings():
        big = dr["rect"].width > 60 and dr["rect"].height > 30
        w = dr.get("width") or 0
        for txt, r in words:
            if dr["type"] in ("f", "fs") and not big:
                inter = fitz.Rect(dr["rect"]) & r
                if not inter.is_empty and inter.width > 0.2 and inter.height > 0.2:
                    tm.append((txt, "fill", [round(v, 1) for v in dr["rect"]]))
                    continue
            if dr["type"] in ("s", "fs") and w:
                for it in dr["items"]:
                    if it[0] == "l" and seg_dist(it[1], it[2], r) < w / 2 - 0.05:
                        tm.append((txt, "line", [round(it[1].x, 1), round(it[1].y, 1), round(it[2].x, 1),
                                                 round(it[2].y, 1)]))
                        break
    d.close()
    return dict(pdf_text_text_overlaps=tt, pdf_text_mark_overlaps=tm, pdf_words=len(words))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE))
    out = pathlib.Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    for k, (p, _) in SOURCES.items():
        assert p is None or p.exists(), p
    L, gates = collect()
    fig, axes, titles = draw(L)
    tg = text_gates(fig)
    pdf, png = out / "Fig3_search_space.pdf", out / "Fig3_search_space.png"
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=200)
    plt.close(fig)
    pg = pdf_gates(pdf)
    po = pdf_overlaps(pdf)
    unused = [k for k in L.index if k not in L.used]
    gates.update(tg)
    gates.update(po)
    gates.update(pdf=pg, ledger_rows=len(L.rows), ledger_values_not_drawn=unused,
                 swarm_worst_overlap=axes["swarm_overlap"])
    sd = out / "Source_Data_Fig3_search_space.csv"
    with open(sd, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["figure", "panel", "row", "field", "value", "source_table",
                                           "specification_label"])
        w.writeheader()
        w.writerows(L.rows)
    checks = {
        "width_exactly_518.4_pt": abs(pg["page_pt"][0] - 518.4) < 0.01,
        "height_le_432_pt": pg["page_pt"][1] <= 432.0,
        "fonts_arial_only": all("Arial" in f for f in pg["fonts"]),
        "min_text_ge_5.7_pt": pg["min_span_size_pt"] >= 5.7 and min(tg["font_sizes"]) >= 5.7,
        "no_text_collisions": not tg["text_collisions"],
        "no_text_text_overlaps_in_pdf": not po["pdf_text_text_overlaps"],
        "no_text_over_marks_in_pdf": not po["pdf_text_mark_overlaps"],
        "titles_7pt_letters_8pt_bold": sorted({round(t.get_fontsize(), 2) for t in titles}) == [FS],
        "no_text_outside_canvas": not tg["text_outside_canvas"],
        "every_ledger_value_drawn_or_printed": not unused,
        "panel_b_swarm_markers_do_not_overlap": max(axes["swarm_overlap"].values()) <= 0.0,
    }
    gates["checks"] = checks
    gates["inputs_sha256"] = {SOURCES[k][1]: sha(p) for k, (p, _) in SOURCES.items() if p is not None}
    gates["inputs_sha256"]["MCP/source_data/Source_Data_Fig3_search_space.csv (read-only cross-check)"] = sha(CURRENT_SD)
    gates["outputs_sha256"] = {x.name: sha(x) for x in (pdf, png, sd)}
    gates["script_sha256"] = sha(__file__)
    print(json.dumps(gates, indent=1, ensure_ascii=False, default=str))
    bad = [k for k, v in checks.items() if not v]
    if bad:
        sys.exit("GATES FAILED: " + ", ".join(bad))
    print("all gates pass")


if __name__ == "__main__":
    main()
