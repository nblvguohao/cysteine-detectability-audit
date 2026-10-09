"""Figure 1 (new): the framework of the paper - where each artefact enters a bottom-up workflow (a) and how the study
is built (b). A schematic, not a data plot.

WHY: the manuscript had no overview figure; both published NC Articles read for style on 2026-09-21 open with one, and
the journal says "schemes are not used; ... experimental procedures should be submitted as figures".

RULES
  R1  Every number drawn is a number the manuscript text (R22) already states; gate G1 checks each against R22 with
      whitespace normalised. No number is introduced here.
  R2  Style from scripts/nc_figure_style_2026-09-21.py: 183 mm wide, drawn at print size, Arial 7 pt, 8 pt bold panel
      letters, every line >= 1 pt (gate G2), only two font sizes (gate G3).
  R3  Text collisions are measured with the tree's detector (gate G4, positive control included).
  R4  No text may fall outside the canvas (gate G5). Added after the second draft lost its last line off the bottom edge,
      which the collision detector cannot see.
INTERPRETER: /path/to/venv2/bin/python
"""
import hashlib, importlib.util, json, pathlib, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "figures_r22"
SELF = pathlib.Path(__file__).read_bytes()


def load(name):
    p = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name[:-3].replace("-", "_"), p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m, hashlib.sha256(p.read_bytes()).hexdigest()


S, S_SHA = load("nc_figure_style_2026-09-21.py")
S.apply()
# Keep text editable in the exported vector files even when this script is run
# independently of the shared style module.
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})

STAGES = ["Sample and\nenrichment", "Protease\ndigestion", "LC-MS/MS\nacquisition", "Database\nsearch",
          "Site table", "Site-level\ncomparison"]
# artefact -> (stage index, short name, what it does)
ARTEFACTS = [
    (4, 0, "Positive = identified", "specific enrichment\nleaves no identified-\nbut-unmodified class"),
    (1, 1, "Cleavage geometry", "K/R spacing mimics\na basic-residue\nmotif"),
    (3, 2, "Protein abundance", "detection tracks\nabundance, not\nmodification"),
    (5, 3, "Search space", "a chemistry the\nsearch never offered\ncannot be called"),
    (2, 4, "Multi-Cys peptides", "inferred sites flip\nneighbour-cysteine\nsigns"),
]


def box(ax, x, y, w, h, text, fc="white", ec=S.INK, lw=S.LW, fs=S.FS, weight="normal", color=S.INK, ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=lw))
    tx = x + w / 2 if ha == "center" else x + 2
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, fontweight=weight, color=color, linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, color=S.INK):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=S.LW, color=color,
                                 shrinkA=0, shrinkB=0))


def badge(ax, x, y, n, r=2.6):
    ax.add_patch(Circle((x, y), r, fc=S.ARTEFACT[n], ec="none"))
    ax.text(x, y, str(n), ha="center", va="center", fontsize=S.FS, fontweight="bold", color="white")


def main():
    width_mm = 183                       # Nature Communications double-column width
    H = 116
    W = width_mm                         # one data unit = 1 mm
    fig = plt.figure(figsize=(width_mm * S.MM, H * S.MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

    # ---------------- a: workflow -------------------------------------------------------------------------------
    ax.text(1, H - 1, "a", fontsize=S.FS_PANEL, fontweight="bold", va="top")
    ax.text(6, H - 1, "Where each artefact enters a bottom-up site-mapping workflow", fontsize=S.FS, va="top")
    n = len(STAGES); bw, gap, y0, bh = 25.5, 4.9, H - 22, 11
    xs = [4 + i * (bw + gap) for i in range(n)]
    for i, (x, name) in enumerate(zip(xs, STAGES)):
        final = i == n - 1
        box(ax, x, y0, bw, bh, name, fc="#f2f2f2" if not final else "white", ec=S.INK)
        if i < n - 1:
            arrow(ax, x + bw + 0.4, y0 + bh / 2, x + bw + gap - 0.4, y0 + bh / 2)
    for num, stage, name, what in ARTEFACTS:
        cx = xs[stage] + bw / 2
        ax.plot([cx, cx], [y0 - 0.2, y0 - 5.2], color=S.ARTEFACT[num], lw=S.LW)
        badge(ax, cx, y0 - 8.2, num)
        ax.text(cx, y0 - 12.4, name, ha="center", va="top", fontsize=S.FS, fontweight="bold", color=S.ARTEFACT[num])
        ax.text(cx, y0 - 12.4 - 3.9 * (name.count("\n") + 1), what, ha="center", va="top", fontsize=S.FS, color=S.INK, linespacing=1.25)
    cx = xs[-1] + bw / 2
    ax.text(cx, y0 - 5.0, "an unmatched\nbackground carries\nall five into the\nresult", ha="center", va="top",
            fontsize=S.FS, color=S.GREY, style="italic", linespacing=1.25)

    # ---------------- b: study design ---------------------------------------------------------------------------
    top = 60
    
    ax.text(1, top + 2, "b", fontsize=S.FS_PANEL, fontweight="bold", va="top")
    ax.text(6, top + 2, "How the study is built", fontsize=S.FS, va="top")

    c1x, c1w = 4, 62
    heads = [("Literature survey", "942 screened  >  367 full text\n> 120 read  >  74 classified"),
             ("Published claims", "51 extracted  >  37 re-testable\n> 35 re-tested  >  28 with a baseline\n(15 papers, 20 site tables)"),
             ("Public data", "9 site-level datasets; deposits\nPXD048216, PXD015307,\nPXD072089 and PXD072035")]
    ys = [top - 14, top - 32, top - 52]
    hs = [12.5, 16.5, 16.5]
    for (t, body), y, h in zip(heads, ys, hs):
        ax.add_patch(FancyBboxPatch((c1x, y), c1w, h, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec=S.INK, lw=S.LW))
        ax.text(c1x + 2.5, y + h - 2.2, t, fontsize=S.FS, fontweight="bold", va="top")
        ax.text(c1x + 2.5, y + h - 6.2, body, fontsize=S.FS, va="top", linespacing=1.3)

    c2x, c2w = 80, 50
    ax.text(c2x, top - 3.5, "One control per artefact", fontsize=S.FS, fontweight="bold", va="top")
    controls = [(1, "background matched on theoretical\ndetectability"), (2, "both sets restricted to\nsingle-cysteine peptides"),
                (3, "abundance-matched protein\nbackground"), (4, "fraction of identified cysteines\nthat carry a site"),
                (5, "both candidates in the\nsearch space")]
    for k, (num, txt) in enumerate(controls):
        y = top - 12 - k * 9.8
        badge(ax, c2x + 2.6, y, num, r=2.4)
        ax.text(c2x + 7, y, txt, fontsize=S.FS, va="center", linespacing=1.2)
    ax.text(c2x, top - 57, "each paired with a same-size random control", fontsize=S.FS, va="top", color=S.GREY, style="italic")

    c3x, c3w = 138, 41
    outs = [("Pre-analysis checklist", "nine items (Table 1)"), ("Label semantics", "three tiers, T1 to T3"),
            ("Benchmark specification", "and what it cannot answer"),
            ("Self-audit", "our ranker on two public\ncohorts, PXD044043 and\nPXD072089")]
    oy = [top - 13, top - 24, top - 35, top - 52]
    oh = [9, 9, 9, 15]
    for (t, b), y, h in zip(outs, oy, oh):
        ax.add_patch(FancyBboxPatch((c3x, y), c3w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    fc="#f2f2f2" if t != "Self-audit" else "white", ec=S.INK, lw=S.LW))
        ax.text(c3x + 2.5, y + h - 2.0, t, fontsize=S.FS, fontweight="bold", va="top")
        ax.text(c3x + 2.5, y + h - 5.6, b, fontsize=S.FS, va="top", linespacing=1.25)
    for y in (ys[0] + hs[0] / 2, ys[1] + hs[1] / 2, ys[2] + hs[2] / 2):
        arrow(ax, c1x + c1w + 0.6, y, c2x - 1.6, y)
    for y in (oy[0] + 4.5, oy[1] + 4.5, oy[2] + 4.5, oy[3] + 7.5):
        arrow(ax, c2x + c2w + 1.0, y, c3x - 0.8, y)

    # ---------------- gates ---------------------------------------------------------------------------------------
    OUT.mkdir(exist_ok=True)
    text_blob = " ".join(t.get_text() for t in fig.findobj(matplotlib.text.Text))
    ms = re.sub(r"\s+", " ", (ROOT / "reports/MANUSCRIPT_SUBMISSION_R22_2026-09-21.md").read_text(encoding="utf-8"))
    nums = sorted(set(re.findall(r"(?<![\w.])\d[\d,]*(?:\.\d+)?", text_blob)) - {"1", "2", "3", "4", "5"})
    missing = [x for x in nums if not re.search(r"(?<![\w.,])" + re.escape(x) + r"(?![\d])", ms)]
    chk, chk_sha = load("check_figure_text_overlaps_2026-09-19.py")
    res = chk.audit_figure("fig1_overview", fig)
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); fw, fh = fig.canvas.get_width_height()
    outside = [t.get_text()[:40] for t in list(ax.texts) + list(fig.texts) if (t.get_text() or "").strip()
               and (lambda e: e.x0 < -0.5 or e.y0 < -0.5 or e.x1 > fw + 0.5 or e.y1 > fh + 0.5)(t.get_window_extent(rend))]
    ctl = plt.figure(); ctl.text(.5, .5, "ctl a", fontsize=10); ctl.text(.5, .5, "ctl b", fontsize=10)
    pc = chk.audit_figure("ctl", ctl)["n_collisions"]; plt.close(ctl)
    gates = dict(G1_numbers_not_in_R22=missing, G1_numbers_checked=nums,
                 G2_min_linewidth=S.min_linewidth(fig), G3_font_sizes=S.font_sizes(fig),
                 G4_collisions=res["n_collisions"], G4_collision_pairs=[(c["text_a"], c["text_b"]) for c in res["collisions"]],
                 G4_positive_control=pc, G5_text_outside_canvas=outside)
    fig.savefig(OUT / "Fig1_overview.pdf"); fig.savefig(OUT / "Fig1_overview.svg")
    fig.savefig(OUT / "Fig1_overview.png", dpi=600); fig.savefig(OUT / "Fig1_overview.tiff", dpi=600)
    plt.close(fig)
    ok = not missing and gates["G2_min_linewidth"] >= 1.0 and set(gates["G3_font_sizes"]) <= {7.0, 8.0} \
        and res["n_collisions"] == 0 and pc > 0 and not outside
    (ROOT / "results/fig1_overview_2026-09-21_audit.json").write_text(json.dumps(dict(
        script=pathlib.Path(__file__).name, script_sha256=hashlib.sha256(SELF).hexdigest(), style_sha256=S_SHA,
        detector_sha256=chk_sha, size_mm=[W, H], gates=gates, all_pass=ok), indent=2, ensure_ascii=False))
    print(json.dumps(gates, indent=1, ensure_ascii=False)); print("all_pass", ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
